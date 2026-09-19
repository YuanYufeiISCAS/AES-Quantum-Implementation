#include "common_model.hpp"
#include <lsqecc/scheduler/wave_scheduler.hpp>
#include <lsqecc/layout/ascii_layout_spec.hpp>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>

namespace lsqecc {
using json = nlohmann::json;
namespace {
constexpr uint32_t TOKEN_BASE = 0x80000000u;
constexpr size_t UNSTARTED = std::numeric_limits<size_t>::max();
const std::map<std::string,size_t> DURATIONS{
    {"cx",2},{"h",3},{"mx",1},{"mz",1},{"reset",1},
    {"x",0},{"z",0},{"ccz_supply",0},{"bell",3},{"z_cond",0},{"x_cond",0},{"barrier",0}};
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
Cell cell(const json& j) { return {j.at(0).get<int>(),j.at(1).get<int>()}; }
json coord(Cell c) { return json::array({c.row,c.col}); }
struct Reservation { size_t end; std::optional<DensePatch> previous; };
struct Stage { std::string kind; size_t start=0,end=0; };
struct Context {
    const json& spec;
    std::ostream& output;
    size_t now=0, emitted=0, makespan=0;
    std::map<Cell,Reservation> reservations;
    std::vector<Cell> locations;
    std::vector<std::string> state, role;
    std::vector<size_t> ends;
    std::vector<std::vector<size_t>> dependencies;
    std::map<size_t,Stage> stages;
    std::map<std::string,size_t> counts;
    std::priority_queue<std::pair<size_t,size_t>,std::vector<std::pair<size_t,size_t>>,
                        std::greater<std::pair<size_t,size_t>>> finishes;
    void complete(size_t index) {
        const auto& op=spec["operations"][index];
        const std::string k=op["kind"];
        if (k=="mz" || k=="mx") state.at(op["qubits"][0])="measured";
        if (k=="reset") state.at(op["qubits"][0])="clean";
        if (k=="ccz_supply") for (auto q:op["qubits"]) state.at(q)="magic";
        if (k=="bell") for (auto q:op["qubits"]) state.at(q)="measured";
        if (k=="cx" && role.at(op["qubits"][1])=="injection_port")
            state.at(op["qubits"][1])="entangled_port";
    }
    void advance(size_t tick, DenseSlice& slice) {
        now=tick;
        for (auto it=reservations.begin();it!=reservations.end();) {
            if (it->second.end>now) { ++it; continue; }
            if (it->second.previous) slice.place_dense_patch_at(it->first,*it->second.previous);
            else slice.clear_patch_at(it->first);
            it=reservations.erase(it);
        }
        while (!finishes.empty() && finishes.top().first<=now) {
            complete(finishes.top().second);
            finishes.pop();
        }
    }
};
thread_local Context* active=nullptr;
struct Scope {
    explicit Scope(Context& ctx) { require(!active,"nested compile");active=&ctx; }
    ~Scope() {active=nullptr;}
};
class Stream final:public LSInstructionStream {
public:
    std::vector<LSInstruction> instructions;
    tsl::ordered_set<PatchId> ids;
    size_t cursor=0;
    LSInstruction get_next_instruction() override {return std::move(instructions.at(cursor++));}
    bool has_next_instruction() const override {return cursor<instructions.size();}
    const tsl::ordered_set<PatchId>& core_qubits() const override {return ids;}
};
InstructionApplicationResult blocked(const std::string& why) {
    return {std::make_unique<std::runtime_error>(why),{}};
}
}

std::ostream& operator<<(std::ostream& out,const CommonModelInstruction& op) {
    return out<<"CommonModelInstruction "<<op.index;
}

InstructionApplicationResult apply_common_model_instruction(
        DenseSlice& slice,LSInstruction& instruction,Router& router) {
    require(active!=nullptr,"missing common-model context");
    auto& c=*active;
    const size_t index=std::get<CommonModelInstruction>(instruction.operation).index;
    if (c.ends[index]!=UNSTARTED) {
        if (c.now<c.ends[index]) return {nullptr,{instruction}};
        return {nullptr,{}};
    }
    const auto& op=c.spec["operations"][index];
    const std::string kind=op["kind"];
    const auto qubits=op["qubits"].get<std::vector<PatchId>>();
    const size_t duration=DURATIONS.at(kind), module=op["module"];
    if (op.contains("guard")) {
        const size_t producer=op["guard"]["bit"];
        if (c.ends[producer]==UNSTARTED || c.now<c.ends[producer]+c.spec.value("feedback_latency",size_t{0}))
            return blocked("measurement feedback not ready");
    }
    if (kind=="z_cond" || kind=="x_cond") for (const auto& bit:op["controls"]) {
        const size_t producer=bit["event"];
        if (c.ends[producer]==UNSTARTED || c.now<c.ends[producer]+c.spec.value("feedback_latency",size_t{0}))
            return blocked("Pauli correction feedback not ready");
    }
    auto& stage=c.stages[module];
    if (duration && c.spec.value("scheduling_policy",std::string("module_homogeneous"))=="module_homogeneous" &&
        stage.end>c.now && (stage.kind!=kind || stage.start!=c.now))
        return blocked("module homogeneous stage active");
    std::vector<Cell> endpoints,footprint;
    if (kind!="barrier") {
        for (auto q:qubits) {
            auto p=c.locations.at(q);
            if (c.reservations.contains(p)) return blocked("operand occupied");
            require(c.state[q]!="measured" || kind=="reset" ||
                    (kind=="ccz_supply" && c.spec.value("resource_reuse",std::string("explicit_reset"))=="accepted_replacement"),
                    "measured operand reused before reset");
            endpoints.push_back(p);
        }
        footprint=endpoints;
    }
    if (kind=="ccz_supply") for (auto q:qubits)
        require(c.role.at(q)=="CCZ_resource" && (c.state.at(q)=="clean" ||
                (c.state.at(q)=="measured" && c.spec.value("resource_reuse",std::string("explicit_reset"))=="accepted_replacement")),
                "CCZ supply must target three clean declared resource sites");
    if (kind=="bell") {
        require(c.role.at(qubits[0])=="injection_port" && c.role.at(qubits[1])=="CCZ_resource",
                "Bell operands must be port, resource");
        require(c.state.at(qubits[0])=="entangled_port" && c.state.at(qubits[1])=="magic",
                "Bell input lifecycle");
        const auto a=endpoints[0], b=endpoints[1];
        require(a.col==b.col && std::abs(a.row-b.row)==1, "Bell requires adjacent vertical pair");
        require(op.at("parity")=="ZZ", "fixed-orientation Bell parity");
    }
    if (kind=="reset") require(c.state.at(qubits[0])=="measured","reset without preceding measurement");
    if (kind=="cx") {
        auto region=router.find_routing_ancilla(slice,qubits[0],PauliOperator::Z,qubits[1],PauliOperator::X);
        if (!region) return blocked("no vertex-disjoint directed route");
        std::vector<Cell> interior;
        for (const auto& p:region->cells) interior.push_back(p.cell);
        std::reverse(interior.begin(),interior.end());
        footprint={endpoints[0]};
        footprint.insert(footprint.end(),interior.begin(),interior.end());
        footprint.push_back(endpoints[1]);
    } else if (kind=="h") {
        const Cell p=endpoints[0];
        const auto corner=c.spec["layout"]["patches"][qubits[0]].value("h_corner",json::array({1,1}));
        const Cell delta=cell(corner);
        require(std::abs(delta.row)==1 && std::abs(delta.col)==1,"invalid H corner");
        footprint={p,{p.row+delta.row,p.col},{p.row,p.col+delta.col},
                     {p.row+delta.row,p.col+delta.col}};
        const auto bound=slice.get_layout().furthest_cell();
        for (size_t j=1;j<footprint.size();++j) {
            const auto f=footprint[j];
            require(f.row>=0 && f.col>=0 && f.row<=bound.row && f.col<=bound.col,"H out of bounds");
            if (!slice.is_cell_free(f)) return blocked("full restoring-H footprint unavailable");
        }
    }
    std::set<Cell> unique(footprint.begin(),footprint.end());
    require(unique.size()==footprint.size(),"duplicate footprint");
    for (auto f:footprint) {
        if (c.reservations.contains(f)) return blocked("footprint occupied");
        if (std::find(endpoints.begin(),endpoints.end(),f)==endpoints.end())
            require(slice.is_cell_free(f),"route traverses live patch");
    }
    json event=op;
    event["start"]=c.now;event["end"]=c.now+duration;
    event["duration"]=duration;event["dependencies"]=c.dependencies[index];
    event["footprint"]=json::array();
    for (auto f:footprint) event["footprint"].push_back(coord(f));
    if (op.contains("guard")) {
        event["feedback_ready"]=c.ends[op["guard"]["bit"].get<size_t>()]+c.spec.value("feedback_latency",size_t{0});
        event["branch_policy"]="reserve_both_outcomes_execute_only_if_predicate";
    }
    c.output<<event.dump()<<'\n';
    c.ends[index]=c.now+duration;
    c.makespan=std::max(c.makespan,c.ends[index]);
    ++c.emitted;++c.counts[kind];
    if (!duration) {c.complete(index);return {nullptr,{}};}
    if (stage.end<=c.now) stage={kind,c.now,c.now+duration};
    for (auto f:footprint) {
        c.reservations.emplace(f,Reservation{c.now+duration,slice.patch_at(f)});
        if (!slice.patch_at(f)) {
            auto patch=DensePatch::from_sparse_patch(LayoutHelpers::basic_square_patch(f));
            patch.type=PatchType::Routing;
            slice.place_dense_patch_at(f,patch);
        }
        slice.set_patch_activity(f,PatchActivity::Reserved);
    }
    c.finishes.emplace(c.now+duration,index);
    return {nullptr,{instruction}};
}

json compile_ccz_model(const json& spec,std::ostream& output) {
    const auto began=std::chrono::steady_clock::now();
    require(spec.at("schema")=="lsc-ccz-template-v1","unknown schema");
    const std::string policy=spec.value("scheduling_policy",std::string("module_homogeneous"));
    require(policy=="module_homogeneous" || policy=="dependency_overlap","unknown scheduling policy");
    const auto& l=spec.at("layout");
    const int rows=l.at("rows"),cols=l.at("cols");
    require(rows>0 && cols>0 && int64_t(rows)*cols<=1000000,"invalid layout");
    std::vector<std::string> grid(rows,std::string(cols,'r'));
    Context ctx{spec,output};
    std::map<Cell,PatchId> located;
    for (const auto& p:l.at("patches")) {
        const PatchId q=p.at("id");const Cell f=cell(p.at("cell"));
        require(q==ctx.locations.size() && q<TOKEN_BASE,"patch IDs must be dense and ordered");
        require(f.row>=0 && f.row<rows && f.col>=0 && f.col<cols,"patch out of bounds");
        require(located.emplace(f,q).second,"duplicate patch location");
        grid[f.row][f.col]='Q';ctx.locations.push_back(f);
        ctx.role.push_back(p.value("role",std::string("data")));
        ctx.state.push_back((ctx.role.back()=="injection_port" || ctx.role.back()=="CCZ_resource")?"clean":"live");
    }
    std::string ascii;
    for (const auto& row:grid) ascii+=row+'\n';
    LayoutFromSpec layout(ascii,DistillationOptions{});
    Stream stream;
    for (const auto& [f,q]:located) stream.ids.insert(q);
    DenseSlice slice(layout,stream.ids);
    Scope scope(ctx);
    const auto& ops=spec.at("operations");
    require(ops.is_array() && ops.size()<=5000000,"invalid operation count");
    ctx.ends.assign(ops.size(),UNSTARTED);ctx.dependencies.resize(ops.size());
    std::vector<size_t> last(ctx.locations.size(),UNSTARTED);
    uint32_t token=TOKEN_BASE;
    for (size_t i=0;i<ops.size();++i) {
        const auto& op=ops[i];const std::string k=op.at("kind");
        require(op.at("id")==i && DURATIONS.contains(k),"invalid operation ID/kind");
        require(!op.contains("duration"),"duration cannot be overridden");
        const auto q=op.at("qubits").get<std::vector<PatchId>>();
        require(k=="barrier" || q.size()==(k=="ccz_supply"?3u:(k=="cx" || k=="bell")?2u:1u),"wrong arity");
        require(std::set<PatchId>(q.begin(),q.end()).size()==q.size(),"duplicate operand");
        LSInstruction instruction{.operation=CommonModelInstruction{i,q},.clients={}};
        std::set<size_t> deps,explicit_deps;
        for (auto wire:q) {
            require(wire<last.size(),"undeclared patch");
            if (last[wire]!=UNSTARTED) deps.insert(last[wire]);
            last[wire]=i;
        }
        for (const auto& dep:op.value("after",json::array())) explicit_deps.insert(dep.get<size_t>());
        if (k=="mx" || k=="mz" || k=="bell") require(op.at("result")==i,"measurement result ID mismatch");
        else require(!op.contains("result"),"non-measurement result");
        if (op.contains("guard")) {
            const size_t producer=op["guard"]["bit"];
            const int value=op["guard"]["equals"];
            require((k=="cx" || k=="z") && (value==0 || value==1),"unsupported predicate");
            require(producer<i && (ops[producer]["kind"]=="mx" || ops[producer]["kind"]=="mz" || ops[producer]["kind"]=="bell"),"invalid result producer");
            explicit_deps.insert(producer);
        }
        if (k=="z_cond" || k=="x_cond") {
            require(op.at("controls").size()==3, "Pauli correction polynomial arity");
            for (const auto& bit:op["controls"]) {
                const size_t producer=bit["event"];
                require(producer<i && ops[producer]["kind"]=="bell", "invalid Pauli result producer");
                explicit_deps.insert(producer);
            }
        }
        for (auto predecessor:explicit_deps) {
            require(predecessor<i && token<UINT32_MAX,"invalid explicit edge");
            stream.instructions[predecessor].clients.insert(token);
            instruction.clients.insert(token++);deps.insert(predecessor);
        }
        ctx.dependencies[i]={deps.begin(),deps.end()};
        stream.instructions.push_back(std::move(instruction));
    }
    for (const auto& instruction:stream.instructions)
        require(instruction.get_patch_dependencies().size()<256,"upstream dependency width exceeded");
    auto router=std::make_unique<CustomDPRouter>();
    const std::string search=spec.value("router",std::string("dijkstra"));
    require(search=="dijkstra" || search=="astar","invalid router");
    router->set_graph_search_provider(search=="astar"?GraphSearchProvider::AStar:GraphSearchProvider::Djikstra);
    DenseSlicingOptions options;options.pipeline_mode=PipelineMode::Wave;
    WaveScheduler scheduler(std::move(stream),layout,std::move(router),options);
    DensePatchComputationResult statistics;
    size_t tick=0, stalled=0;
    while (!scheduler.done()) {
        require(tick<=spec.value("max_cycles",size_t{1000000}),"logical cycle budget exceeded");
        require(std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<
                spec.value("timeout_seconds",7200.0),"wall clock budget exceeded");
        ctx.advance(tick,slice);
        const size_t before=ctx.emitted;
        auto stats=scheduler.schedule_wave(slice,[](const auto&){},statistics);
        if (ctx.emitted==before && ctx.finishes.empty()) ++stalled;else stalled=0;
        require(stalled<spec.value("max_stalled_cycles",size_t{100}),"deadlock: "+stats.blocked_cause);
        if (tick%1000==0) std::cerr<<"cycle="<<tick<<" completed_starts="<<ctx.emitted<<"/"<<ops.size()<<'\n';
        ++tick;
    }
    ctx.advance(ctx.makespan,slice);
    require(ctx.emitted==ops.size() && ctx.reservations.empty() && ctx.finishes.empty(),"incomplete schedule");
    for (size_t q=0;q<ctx.role.size();++q)
        if (ctx.role[q]=="injection_port") require(ctx.state[q]=="clean","unclean final injection port");
        else if (ctx.role[q]=="CCZ_resource")
            require(ctx.state[q]=="clean" || (ctx.state[q]=="measured" &&
                    spec.value("resource_reuse",std::string("explicit_reset"))=="accepted_replacement"),
                    "unclean final resource site");
    return {{"schema","lsc-ccz-template-result-v1"},{"status","compiled_requires_independent_audit"},
        {"latency_logical_cycles",ctx.makespan},{"reserved_patches",size_t(rows)*cols},
        {"reserved_patch_cycles",size_t(rows)*cols*ctx.makespan},
        {"operations",ctx.emitted},{"primitive_counts",ctx.counts},{"final_patch_states",ctx.state},
        {"scheduling_policy",policy},{"router",search},
        {"compiler_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()}};
}
}
