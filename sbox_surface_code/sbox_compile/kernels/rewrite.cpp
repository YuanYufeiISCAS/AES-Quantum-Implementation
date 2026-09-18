// Whole-matrix native-VDP search with deterministic proposal and route budgets.
// GL(4,2) routines are from the authors' Mixcolumn/compiler_baselines/refine_gl4.h.
// The free-output identity follows refine_free_output.h: F=S M T, C=Q M,
// S'=Q S Q^-1, hence S' C T=Q F. Output permutations are *proposals*, not moves.
#include "gl4.hpp"
#include <cmath>
#include <ctime>
#include <deque>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>

namespace joint {
using Gate = mixcolumn_refine_gl4::Gate;
using Word = mixcolumn_refine_gl4::Word;
using Rows = std::vector<std::uint64_t>;
using Perm = std::vector<int>;
namespace gl4 = mixcolumn_refine_gl4;

struct Request {
    int n=0, rows=0, cols=0, height=0, width=0;
    std::uint64_t seed=0;

    bool free_output=false, has_seed=false;
    int keep=4, max_proposals=0, max_routes=0, max_gates=2048;
    Rows matrix;
    Word seed_word;
    Perm seed_q;
};
struct RoutedGate { Gate gate; std::vector<int> path; };
using Layer = std::vector<RoutedGate>;
struct State {
    Word word;
    Perm q;
    std::vector<Layer> layers;
    std::string branch, move;
    int origin_ordinal=-1;
    Word origin_word;
    Perm origin_q;
    int accepted_walk_steps=0;
};

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
Perm identity(int n) { Perm p(n); std::iota(p.begin(),p.end(),0); return p; }
Perm inverse(const Perm& p) {
    Perm q(p.size(),-1);
    for (int i=0;i<int(p.size());++i) {
        require(p[i]>=0 && p[i]<int(p.size()) && q[p[i]]<0,"bad permutation");
        q[p[i]]=i;
    }
    return q;
}
Rows replay(int n, const Word& word) {
    Rows rows(n);
    for (int i=0;i<n;++i) rows[i]=std::uint64_t(1)<<i;
    for (auto gate: word) {
        require(gate.first>=0 && gate.first<n && gate.second>=0 && gate.second<n
                && gate.first!=gate.second,"bad word endpoints");
        rows[gate.second]^=rows[gate.first];
    }
    return rows;
}
bool realizes(const Request& request, const State& state) {
    (void)inverse(state.q);
    if (!request.free_output && state.q!=identity(request.n)) return false;
    const auto actual=replay(request.n,state.word);
    for (int i=0;i<request.n;++i) if (actual[i]!=request.matrix[state.q[i]]) return false;
    return true;
}
bool commute(Gate a, Gate b) { return a.first!=b.second && a.second!=b.first; }
bool cancel(Word& word) {
    bool changed=false;
    for (std::size_t i=0;i<word.size();) {
        bool erased=false;
        for (std::size_t j=i+1;j<word.size();++j) {
            if (word[j]==word[i]) {
                word.erase(word.begin()+j); word.erase(word.begin()+i);
                changed=erased=true; break;
            }
            if (!commute(word[i],word[j])) break;
        }
        if (!erased) ++i;
    }
    return changed;
}
std::uint64_t key(const State& state, bool multiset=false) {
    std::uint64_t hash=1469598103934665603ULL;
    auto add=[&](std::uint64_t v) { hash^=v; hash*=1099511628211ULL; };
    for (int v:state.q) add(v+1);
    add(257);
    Word word=state.word;
    if (multiset) std::sort(word.begin(),word.end());
    for (auto gate:word) { add(gate.first+1); add(gate.second+1); }
    return hash;
}
auto score(const State& state) {
    return std::make_tuple(2*int(state.layers.size()),int(state.word.size()),key(state));
}

// Genuine synthesis: row additions eliminate the requested matrix itself.
// Fixed-output mode selects columns in randomized order and uses one additive
// pivot instead of a three-CNOT row swap. The transposed variant is converted
// back by reversing the word and swapping every directed CNOT endpoint.
State fresh(const Request& request, std::mt19937_64& rng, int ordinal) {
    State result;
    result.branch="fresh_whole_matrix_gaussian";
    result.origin_ordinal=ordinal;
    result.q=identity(request.n);
    Rows reduced=request.matrix;
    const bool free=request.free_output && ordinal%4==1;
    const bool transpose=!free && ordinal%2==1;
    if (transpose) {
        Rows transposed(request.n,0);
        for (int i=0;i<request.n;++i)
            for(int j=0;j<request.n;++j) transposed[j]|=((reduced[i]>>j)&1ULL)<<i;
        reduced=std::move(transposed);
    }
    Perm columns=identity(request.n), completed(request.n,0);
    if (ordinal) std::shuffle(columns.begin(),columns.end(),rng);
    Word elimination;
    auto add=[&](int c,int t) { reduced[t]^=reduced[c]; elimination.emplace_back(c,t); };
    // Small, recorded algebraic preconditioning opens different factorizations.
    if (ordinal>=4 && ordinal%3==0) {
        for (int k=0;k<1+ordinal%3;++k) {
            int c=int(rng()%request.n),t=int(rng()%request.n);
            if(c!=t) add(c,t);
        }
    }
    for (int column:columns) {
        std::vector<int> pivots;
        for(int r=0;r<request.n;++r)
            if(!completed[r] && ((reduced[r]>>column)&1)) pivots.push_back(r);
        require(!pivots.empty(),"singular fresh matrix");
        int pivot;
        if (free) {
            pivot=pivots[std::size_t(rng()%pivots.size())];
        } else {
            pivot=column;
            if(!((reduced[pivot]>>column)&1)) {
                int source=pivots[std::size_t(rng()%pivots.size())];
                if (ordinal%3==0) {
                    source=*std::min_element(pivots.begin(),pivots.end(),[&](int a,int b) {
                        return __builtin_popcountll(reduced[pivot]^reduced[a])
                               < __builtin_popcountll(reduced[pivot]^reduced[b]);
                    });
                }
                add(source,pivot);
            }
        }
        Perm targets=identity(request.n);
        if(ordinal) std::shuffle(targets.begin(),targets.end(),rng);
        for(int row:targets) if(row!=pivot && ((reduced[row]>>column)&1)) add(pivot,row);
        completed[pivot]=1;
    }
    result.word=Word(elimination.rbegin(),elimination.rend());
    if (free) {
        Perm p(request.n,-1);
        for(int row=0;row<request.n;++row) {
            require(__builtin_popcountll(reduced[row])==1,"free elimination is not a permutation");
            p[row]=__builtin_ctzll(reduced[row]);
        }
        // L A=R. R^-1 L^-1 R=(R^-1) A, with endpoints mapped by p.
        for(auto& gate:result.word) { gate.first=p[gate.first]; gate.second=p[gate.second]; }
        result.q=inverse(p);
    } else {
        for(int row=0;row<request.n;++row)
            require(reduced[row]==(std::uint64_t(1)<<row),"fixed elimination did not reach I");
    }
    if (transpose) {
        std::reverse(result.word.begin(),result.word.end());
        for(auto& gate:result.word) std::swap(gate.first,gate.second);
    }
    cancel(result.word);
    require(realizes(request,result),"fresh whole matrix replay failed");
    result.origin_word=result.word; result.origin_q=result.q;
    result.move=free?"fresh_free_output_gaussian":transpose?"fresh_transposed_gaussian":"fresh_direct_gaussian";
    return result;
}

class Router {
    const Request& request_;
    std::vector<int> occupied_,positions_;
public:
    explicit Router(const Request& request):request_(request),
        occupied_(request.height*request.width,0),positions_(request.n) {
        for(int i=0;i<request.n;++i) {
            int p=(2*(i/request.cols)+1)*request.width+2*(i%request.cols)+1;
            positions_[i]=p; occupied_[p]=1;
        }
    }
    std::vector<int> path(Gate gate,const std::vector<int>& used,std::uint64_t salt) const {
        const int start=positions_[gate.first],goal=positions_[gate.second];
        if(used[start] || used[goal]) return {};
        std::vector<int> parent(occupied_.size(),-1),queue;
        parent[start]=start; queue.push_back(start);
        const std::array<std::pair<int,int>,4> offsets{{{-1,0},{0,-1},{0,1},{1,0}}};
        for(std::size_t head=0;head<queue.size();++head) {
            const int current=queue[head],r=current/request_.width,c=current%request_.width;
            if(current==goal) {
                std::vector<int> result;
                for(int p=goal;;p=parent[p]) { result.push_back(p); if(p==start) break; }
                std::reverse(result.begin(),result.end()); return result;
            }
            std::array<int,4> order{{0,1,2,3}};
            std::rotate(order.begin(),order.begin()+salt%4,order.end());
            if(salt&4) std::reverse(order.begin(),order.end());
            std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
                auto distance=[&](int i) { return std::abs(r+offsets[i].first-goal/request_.width)
                       +std::abs(c+offsets[i].second-goal%request_.width); };
                return distance(a)<distance(b);
            });
            for(int i:order) {
                int nr=r+offsets[i].first,nc=c+offsets[i].second;
                if(nr<0 || nr>=request_.height || nc<0 || nc>=request_.width) continue;
                int next=nr*request_.width+nc;
                if(parent[next]>=0 || used[next] || (occupied_[next] && next!=goal)) continue;
                if(current==start && nc!=c) continue; // control exits vertically
                if(next==goal && r!=goal/request_.width) continue; // target enters horizontally
                parent[next]=current; queue.push_back(next);
            }
        }
        return {};
    }
    bool route(State& state,std::mt19937_64& rng,int policy) const {
        const int count=int(state.word.size());
        std::vector<std::vector<int>> successors(count);
        std::vector<int> indegree(count,0),height(count,1),done(count,0);
        // Exact noncommutation DAG, not the unnecessarily serial same-wire DAG.
        for(int i=count-1;i>=0;--i) for(int j=i+1;j<count;++j)
            if(!commute(state.word[i],state.word[j])) {
                successors[i].push_back(j); ++indegree[j]; height[i]=std::max(height[i],height[j]+1);
            }
        state.layers.clear();
        Word scheduled;
        for(int remaining=count;remaining;) {
            std::vector<int> ready;
            for(int i=0;i<count;++i) if(!done[i] && indegree[i]==0) ready.push_back(i);
            require(!ready.empty(),"CNOT dependency graph contains a cycle");
            if(policy) std::shuffle(ready.begin(),ready.end(),rng);
            auto distance=[&](int i) {
                auto gate=state.word[i]; int a=positions_[gate.first],b=positions_[gate.second];
                return std::abs(a/request_.width-b/request_.width)+std::abs(a%request_.width-b%request_.width);
            };
            std::stable_sort(ready.begin(),ready.end(),[&](int a,int b) {
                if(policy%3==1) return std::make_tuple(-height[a],distance(a))<std::make_tuple(-height[b],distance(b));
                if(policy%3==2) return distance(a)<distance(b);
                return height[a]>height[b];
            });
            std::vector<int> used(occupied_.size(),0),selected;
            Layer layer;
            for(int index:ready) {
                auto route=path(state.word[index],used,rng());
                if(route.empty()) continue;
                for(int p:route) used[p]=1;
                layer.push_back({state.word[index],std::move(route)}); selected.push_back(index);
            }
            if(selected.empty()) return false;
            for(int index:selected) {
                done[index]=1; --remaining; scheduled.push_back(state.word[index]);
                for(int next:successors[index]) --indegree[next];
            }
            state.layers.push_back(std::move(layer));
        }
        require(replay(request_.n,state.word)==replay(request_.n,scheduled),"router changed matrix");
        state.word=std::move(scheduled);
        return true;
    }
};

State replace_window(const State& base,int begin,int length,const Perm& active,
                     const Word& replacement,const gl4::Permutation& q) {
    State result=base;
    Perm endpoint=identity(int(base.q.size()));
    for(int i=0;i<int(active.size());++i) {
        endpoint[active[q[i]]]=active[i]; result.q[active[i]]=base.q[active[q[i]]];
    }
    result.word.assign(base.word.begin(),base.word.begin()+begin);
    for(auto gate:replacement) result.word.emplace_back(active[gate.first],active[gate.second]);
    for(std::size_t i=begin+length;i<base.word.size();++i)
        result.word.emplace_back(endpoint[base.word[i].first],endpoint[base.word[i].second]);
    result.layers.clear();
    return result;
}

void emit_word(const Word& word) {
    std::cout<<'[';
    for(std::size_t i=0;i<word.size();++i) {
        if(i) std::cout<<',';
        std::cout<<'['<<word[i].first<<','<<word[i].second<<']';
    }
    std::cout<<']';
}
void emit_perm(const Perm& p) {
    std::cout<<'['; for(std::size_t i=0;i<p.size();++i) { if(i) std::cout<<','; std::cout<<p[i]; } std::cout<<']';
}
void emit(const Request& request,const State& state) {
    std::cout<<"{\"candidate\":{\"n\":"<<request.n<<",\"layout\":{\"data_rows\":"<<request.rows
        <<",\"data_cols\":"<<request.cols<<",\"grid_rows\":"<<request.height<<",\"grid_cols\":"<<request.width
        <<"},\"mode\":\"vdp\",\"target_rows_hex\":[";
    for(int i=0;i<request.n;++i) { if(i) std::cout<<','; std::cout<<"\"0x"<<std::hex<<request.matrix[i]<<std::dec<<'"'; }
    std::cout<<"],\"output_permutation\":"; emit_perm(state.q);
    std::cout<<",\"layers\":[";
    for(std::size_t i=0;i<state.layers.size();++i) {
        if(i) std::cout<<',';
        std::cout<<"{\"mode\":\"vdp\",\"logical_depth\":2,\"operations\":[";
        for(std::size_t j=0;j<state.layers[i].size();++j) {
            if(j) std::cout<<',';
            const auto& gate=state.layers[i][j];
            std::cout<<"{\"control\":"<<gate.gate.first<<",\"target\":"<<gate.gate.second<<",\"path\":[";
            for(std::size_t k=0;k<gate.path.size();++k) {
                if(k) std::cout<<',';
                std::cout<<'['<<gate.path[k]/request.width<<','<<gate.path[k]%request.width<<']';
            }
            std::cout<<"]}";
        }
        std::cout<<"]}";
    }
    std::cout<<"],\"stats\":{\"surface_depth\":"<<2*state.layers.size()<<",\"layers\":"<<state.layers.size()
        <<",\"cnots\":"<<state.word.size()<<"},\"construction\":\"joint_matrix_native_vdp\"},\"evidence\":{\"branch\":\""
        <<state.branch<<"\",\"origin_ordinal\":"<<state.origin_ordinal<<",\"origin_word\":";
    emit_word(state.origin_word); std::cout<<",\"origin_output_permutation\":"; emit_perm(state.origin_q);
    std::cout<<",\"last_move\":\""<<state.move<<"\",\"accepted_walk_steps\":"<<state.accepted_walk_steps<<"}}";
}

Request read_request() {
    Request request;
    std::string magic;
    std::cin>>magic;
    require(magic=="JOINT_LINEAR_V1","wrong request schema");
    int free=0,has_seed=0;
    std::cin>>request.n>>request.rows>>request.cols>>request.height>>request.width
        >>request.seed>>free>>request.keep>>request.max_proposals>>request.max_routes>>request.max_gates;
    require(request.n>=1 && request.n<=32 && request.rows>=1 && request.rows<=5 && request.cols>=1
        && request.n<=request.rows*request.cols && request.height==2*request.rows+1 && request.width==2*request.cols+1,
        "invalid occupied lattice geometry");
    require((free==0 || free==1) && request.keep>=1 && request.keep<=4
        && request.max_proposals>=1 && request.max_proposals<=64000000
        && request.max_routes>=1 && request.max_routes<=2000000 && request.max_gates>=1 && request.max_gates<=2048,
        "invalid finite search limits");
    request.free_output=free;
    request.matrix.resize(request.n);
    for(auto& row:request.matrix) {
        std::cin>>row; require(row<(std::uint64_t(1)<<request.n),"matrix row outside width");
    }
    std::cin>>has_seed;
    require(has_seed==0 || has_seed==1,"invalid seed presence");
    request.has_seed=has_seed;
    if(request.has_seed) {
        request.seed_q.resize(request.n);
        for(int& value:request.seed_q) std::cin>>value;
        (void)inverse(request.seed_q);
        int count=-1; std::cin>>count;
        require(count>=0 && count<=request.max_gates,"seed word exceeds bounded input");
        request.seed_word.resize(count);
        for(auto& gate:request.seed_word) std::cin>>gate.first>>gate.second;
        State seed; seed.word=request.seed_word; seed.q=request.seed_q;
        require(realizes(request,seed),"seed matrix does not realize request");
    }
    require(bool(std::cin),"truncated request");
    std::string extra; require(!(std::cin>>extra),"trailing request bytes");
    return request;
}

int run(const Request& request) {
    const double start=double(std::clock())/CLOCKS_PER_SEC;
    auto elapsed=[&]() { return double(std::clock())/CLOCKS_PER_SEC-start; };
    std::mt19937_64 rng(request.seed);
    Router router(request);
    std::vector<State> population;
    std::unordered_set<std::uint64_t> seen;
    int proposals=0,routes=0,fresh_count=0,rewrites=0,walks=0;
    int best_depth=std::numeric_limits<int>::max(),best_gates=request.max_gates;
    State walk;
    bool has_walk=false;
    auto expired=[&]() { return proposals>=request.max_proposals || routes>=request.max_routes; };
    auto consider=[&](State candidate,bool force=false) {
        ++proposals;
        if(candidate.word.size()>std::size_t(request.max_gates)) return;
        require(realizes(request,candidate),"proposed algebraic word does not realize Q A");
        if(!force && !seen.insert(key(candidate)).second) return;
        if(seen.size()>50000) seen.clear(); // bounded heuristic tabu, never a proof cache
        if(expired() && !force) return;
        ++routes;
        if(!router.route(candidate,rng,routes%4)) return;
        require(realizes(request,candidate),"routed matrix differs from target");
        const int depth=2*int(candidate.layers.size());
        const bool better=depth<best_depth || (depth==best_depth && int(candidate.word.size())<best_gates);
        if(better) { best_depth=depth; best_gates=int(candidate.word.size()); }
        // A separate exploratory walk can accept equal/longer networks. The
        // retained depth incumbent is never overwritten by an uphill step.
        if(better || !has_walk || (depth<=best_depth+4 && int(candidate.word.size())<=best_gates+4
                                 && (depth<=2*int(walk.layers.size()) || rng()%8==0))) {
            candidate.accepted_walk_steps+=1; walk=candidate; has_walk=true; ++walks;
        }
        const auto duplicate=std::find_if(population.begin(),population.end(),[&](const State& old) {
            return key(old)==key(candidate) && old.q==candidate.q && old.word==candidate.word;
        });
        if(duplicate!=population.end()) {
            if(score(candidate)<score(*duplicate)) *duplicate=std::move(candidate);
        } else population.push_back(std::move(candidate));
        std::sort(population.begin(),population.end(),[](const State& a,const State& b){return score(a)<score(b);});
        if(population.size()>8) {
            // Always keep best; prefer distinct output permutations and gate
            // multisets for the other seven independently feasible members.
            std::vector<State> kept;
            std::set<Perm> qs;
            std::unordered_set<std::uint64_t> sets;
            while(!population.empty() && kept.size()<8) {
                std::size_t selected=0;
                for(std::size_t i=1;i<population.size() && !kept.empty();++i) {
                    auto diversity=[&](std::size_t j) {
                        return std::make_tuple(qs.count(population[j].q)?1:0,sets.count(key(population[j],true))?1:0,score(population[j]));
                    };
                    if(diversity(i)<diversity(selected)) selected=i;
                }
                qs.insert(population[selected].q); sets.insert(key(population[selected],true));
                kept.push_back(std::move(population[selected])); population.erase(population.begin()+selected);
            }
            population=std::move(kept);
        }
    };
    if(request.has_seed) {
        State seed;
        seed.word=request.seed_word; seed.q=request.seed_q;
        seed.origin_word=seed.word; seed.origin_q=seed.q;
        seed.branch="source_bound_network"; seed.move="seed_route_rebuild";
        consider(seed,true);
        State reduced=seed;
        if(cancel(reduced.word)) { reduced.move="commuting_cancellation"; consider(std::move(reduced)); }
    }
    // At least one true matrix-only synthesis is attempted even with a seed.
    consider(fresh(request,rng,fresh_count++),true);
    (void)gl4::table();
    for(int round=0;round<8*request.max_proposals && !expired();++round) {
        if(round%16==0 || population.empty()) {
            consider(fresh(request,rng,fresh_count++));
            if(population.empty()) continue;
        }
        State base=(has_walk && round%3 ? walk:population[rng()%population.size()]);
        if(base.word.empty()) {
            if(request.free_output && fresh_count<16) { consider(fresh(request,rng,fresh_count++)); continue; }
            break;
        }
        if(round%32==31) { has_walk=false; continue; }
        if(round%9==0 && base.word.size()>1) {
            int index=int(rng()%(base.word.size()-1));
            if(commute(base.word[index],base.word[index+1])) {
                std::swap(base.word[index],base.word[index+1]); base.move="adjacent_commuting_walk";
                consider(base);
            }
        }
        const int begin=int(rng()%base.word.size());
        Perm support;
        int length=0;
        for(int end=begin;end<int(base.word.size()) && end<begin+12;++end) {
            Perm next=support;
            for(int wire:{base.word[end].first,base.word[end].second})
                if(std::find(next.begin(),next.end(),wire)==next.end()) next.push_back(wire);
            if(next.size()>4) break;
            support=std::move(next); length=end-begin+1;
        }
        if(!length) continue;
        // Padding wires may participate transiently; the exact local 4x4
        // action is identity on them. They are never removed as obstacles.
        for(int wire=0;support.size()<4 && wire<request.n;++wire)
            if(std::find(support.begin(),support.end(),wire)==support.end()) support.push_back(wire);
        Perm lookup(request.n,-1);
        const int arity=int(support.size());
        for(int i=0;i<arity;++i) lookup[support[i]]=i;
        gl4::Code local=gl4::identity;
        for(int i=begin;i<begin+length;++i)
            local=gl4::append_cnot(local,lookup[base.word[i].first],lookup[base.word[i].second]);
        std::vector<gl4::Permutation> permutations{{0,1,2,3}};
        if(request.free_output) {
            gl4::Permutation p{0,1,2,3};
            do {
                bool physical=true;
                for(int i=arity;i<4;++i) if(p[i]!=i) physical=false;
                if(physical && p!=gl4::Permutation{0,1,2,3}) permutations.push_back(p);
            }
            while(std::next_permutation(p.begin(),p.end()));
            std::shuffle(permutations.begin()+1,permutations.end(),rng);
            permutations.resize(std::min(permutations.size(),std::size_t(1+std::min(5,round%24+1))));
        }
        for(const auto& q:permutations) {
            if(expired()) break;
            const auto target=gl4::left_row_permute(local,q);
            auto variants=gl4::shortest_variants(target,rng(),4);
            if(round%3==0) {
                // One extra trailing generator often exposes a depth-better
                // factorization with the same or slightly larger CNOT count.
                int c=int(rng()%4),t=int(rng()%3); if(t>=c) ++t;
                auto longer=gl4::table().sampled_shortest(gl4::append_cnot(target,c,t),rng);
                longer.emplace_back(c,t); variants.push_back(std::move(longer));
            }
            for(const auto& replacement:variants) {
                if(expired()) break;
                if(int(replacement.size())>length+4) continue;
                bool physical=true;
                for(const auto& gate:replacement) if(gate.first>=arity || gate.second>=arity) physical=false;
                if(!physical) continue; // no virtual fourth patch on a tiny toy
                State candidate=replace_window(base,begin,length,support,replacement,q);
                cancel(candidate.word);
                if(candidate.word==base.word && candidate.q==base.q) continue;
                candidate.move=q==gl4::Permutation{0,1,2,3}?"exact_GL4_rewrite":"exact_freeQ_GL4_suffix_conjugation";
                ++rewrites; consider(std::move(candidate));
            }
        }
    }
    // Selection remains a bounded diversity portfolio, not a single scalar.
    std::sort(population.begin(),population.end(),[](const State& a,const State& b){return score(a)<score(b);});
    std::vector<State> selected;
    std::set<Perm> qs;
    std::unordered_set<std::uint64_t> sets;
    while(!population.empty() && selected.size()<std::size_t(request.keep)) {
        std::size_t index=0;
        for(std::size_t i=1;i<population.size() && !selected.empty();++i) {
            auto value=[&](std::size_t j) {
                return std::make_tuple(qs.count(population[j].q)?1:0,sets.count(key(population[j],true))?1:0,score(population[j]));
            };
            if(value(i)<value(index)) index=i;
        }
        qs.insert(population[index].q); sets.insert(key(population[index],true));
        selected.push_back(std::move(population[index])); population.erase(population.begin()+index);
    }
    std::cout<<"{\"schema\":\"joint-linear-engine-v1\",\"candidates\":[";
    for(std::size_t i=0;i<selected.size();++i) { if(i) std::cout<<','; emit(request,selected[i]); }
    std::cout<<"],\"telemetry\":{\"cpu_seconds\":"<<elapsed()<<",\"proposals\":"<<proposals
        <<",\"routed_words\":"<<routes<<",\"fresh_matrix_starts\":"<<fresh_count<<",\"algebraic_rewrites\":"<<rewrites
        <<",\"walk_admissions\":"<<walks<<",\"threads\":1,\"gl4_states\":"<<gl4::table().state_count()
        <<",\"termination\":\""<<(routes>=request.max_routes?"route_cap_unknown":
            proposals>=request.max_proposals?"proposal_cap_unknown":"finite_search_complete_not_optimal")<<"\"}}\n";
    return 0;
}
} // namespace joint

int main(int argc,char** argv) {
    try {
        joint::require(argc==2 && std::string(argv[1])=="--stdin-v1","usage: joint_linear_backend --stdin-v1");
        return joint::run(joint::read_request());
    } catch(const std::exception& error) {
        std::cerr<<"joint-linear: "<<std::string(error.what()).substr(0,1024)<<'\n'; return 2;
    }
}
