"""Independent streaming space-time audit; does not import the scheduler/router."""
from collections import Counter
import json

DURATIONS={"cx":2,"h":3,"mx":1,"mz":1,"reset":1,"x":0,"z":0,"ccz_supply":0,"bell":3,"z_cond":0,"x_cond":0,"barrier":0}

def require(ok, message):
    if not ok:
        raise ValueError(message)


def expected_dependencies(spec):
    last={}
    result=[]
    for i,op in enumerate(spec["operations"]):
        require(op["id"]==i,"operation index")
        deps=set(op.get("after",[]))
        for q in op["qubits"]:
            if q in last:
                deps.add(last[q])
            last[q]=i
        if "guard" in op:
            deps.add(op["guard"]["bit"])
        if op["kind"] in ("z_cond","x_cond"):
            deps.update(b["event"] for b in op["controls"])
        require(all(0<=j<i for j in deps),"forward dependency")
        result.append(deps)
    return result


def audit(spec, trace_path, result):
    ops=spec["operations"]
    dependencies=expected_dependencies(spec)
    layout=spec["layout"]
    positions={p["id"]:tuple(p["cell"]) for p in layout["patches"]}
    roles={p["id"]:p["role"] for p in layout["patches"]}
    occupied=set(positions.values())
    require(len(occupied)==len(positions),"duplicate static patch")
    state={q:"clean" if roles[q] in ("injection_port","CCZ_resource") else "live" for q in positions}
    ends=[-1]*len(ops)
    starts=[-1]*len(ops)
    reserved={}  # cell -> last positive reservation end; events are start-ordered
    stages={}
    counts=Counter()
    module_times={}
    previous_start=0
    with trace_path.open() as source:
        for line in source:
            event=json.loads(line)
            i=event["id"]
            require(0<=i<len(ops) and ends[i]<0,"duplicate or unknown event")
            op=ops[i]
            for key,value in op.items():
                require(event.get(key)==value,f"changed operation field {i}.{key}")
            kind=op["kind"];q=op["qubits"]
            start,end=event["start"],event["end"]
            require(type(start) is int and type(end) is int and start>=previous_start,"event clock order")
            previous_start=start
            require(event["duration"]==DURATIONS[kind] and end==start+DURATIONS[kind],"primitive duration")
            require(set(event["dependencies"])==dependencies[i],"wrong exported dependencies")
            require(all(0<=ends[j]<=start for j in dependencies[i]),f"dependency violation {i}")
            if "guard" in op:
                bit=op["guard"]["bit"]
                require(kind in ("cx","z") and op["guard"]["equals"] in (0,1),"guard syntax")
                require(ops[bit]["kind"] in ("mx","mz","bell"),"result producer")
                ready=ends[bit]+spec["feedback_latency"]
                require(start>=ready and event["feedback_ready"]==ready,"feedback readiness")
                require(event["branch_policy"]=="reserve_both_outcomes_execute_only_if_predicate","guard erased")
            if kind in ("z_cond","x_cond"):
                require(len(op["controls"])==3,"Pauli polynomial arity")
                for b in op["controls"]:
                    require(ops[b["event"]]["kind"]=="bell","Pauli outcome producer")
                    require(start>=ends[b["event"]]+spec["feedback_latency"],"Pauli feedback")
            footprint=list(map(tuple,event["footprint"]))
            require(len(footprint)==len(set(footprint)),"duplicate footprint")
            endpoints={positions[j] for j in q} if kind!="barrier" else set()
            require(endpoints<=set(footprint),"missing endpoint")
            for r,c in footprint:
                require(0<=r<layout["rows"] and 0<=c<layout["cols"],"footprint out of bounds")
                require(reserved.get((r,c),0)<=start,f"space-time collision at {(r,c)}")
            require(not (set(footprint)-endpoints)&occupied,"route/H crosses live static patch")
            if kind=="cx":
                require(len(footprint)>=3 and footprint[0]==positions[q[0]] and
                        footprint[-1]==positions[q[1]],"CNOT endpoints")
                require(all(abs(a[0]-b[0])+abs(a[1]-b[1])==1
                            for a,b in zip(footprint,footprint[1:])),"discontinuous route")
                require(footprint[1][1]==footprint[0][1] and
                        footprint[-2][0]==footprint[-1][0],"CNOT boundary orientation")
            elif kind=="h":
                r,c=positions[q[0]]
                dr,dc=layout["patches"][q[0]]["h_corner"]
                require(set(footprint)=={(r,c),(r+dr,c),(r,c+dc),(r+dr,c+dc)},"full H footprint")
            else:
                require(set(footprint)==endpoints,"unexpected primitive footprint")
            if kind!="barrier":
                require(all(state[j]!="measured" or kind=="reset" or
                    (kind=="ccz_supply" and spec.get("resource_reuse")=="accepted_replacement") for j in q),"reuse before reset")
            if kind=="ccz_supply":
                require(len(q)==3 and len(set(q))==3,"CCZ triple arity")
                require(all(roles[j]=="CCZ_resource" and (state[j]=="clean" or
                    (state[j]=="measured" and spec.get("resource_reuse")=="accepted_replacement")) for j in q),"illegal CCZ supply")
                for j in q: state[j]="magic"
            elif kind=="bell":
                require(len(q)==2 and roles[q[0]]=="injection_port" and roles[q[1]]=="CCZ_resource","Bell roles")
                a,b=positions[q[0]],positions[q[1]]
                require(a[1]==b[1] and abs(a[0]-b[0])==1 and op["parity"]=="ZZ","Bell adjacency/orientation")
                require(state[q[0]]=="entangled_port" and state[q[1]]=="magic","Bell lifecycle")
                require(op["result"]==i,"Bell result")
                for j in q: state[j]="measured"
            elif kind=="cx" and roles[q[1]]=="injection_port":
                state[q[1]]="entangled_port"
            elif kind in ("mz","mx"):
                require(op["result"]==i,"measurement identifier")
                state[q[0]]="measured"
            elif kind=="reset":
                require(state[q[0]]=="measured","reset without readout")
                state[q[0]]="clean"
            if end>start:
                module=op["module"]
                if spec["scheduling_policy"]=="module_homogeneous":
                    if module in stages and stages[module][1]>start:
                        require(stages[module]==(start,end,kind),"mixed/overlapping module stage")
                    else:
                        stages[module]=(start,end,kind)
                for p in footprint:
                    reserved[p]=end
            ends[i]=end
            starts[i]=start
            counts[kind]+=1
            m=op["module"]
            bounds=module_times.setdefault(m,[start,end])
            bounds[0]=min(bounds[0],start);bounds[1]=max(bounds[1],end)
    require(all(e>=0 for e in ends),"missing events")
    for gadget in spec["gadgets"]:
        if gadget["kind"]=="qand_dagger": continue
        block=ops[gadget["start"]:gadget["end"]]
        supplies=[o for o in block if o["kind"]=="ccz_supply"]
        require(len(supplies)==1,"CCZ supply inventory")
        injection=[o for o in block if o["kind"]=="cx" and roles[o["qubits"][1]]=="injection_port"]
        require(len(injection)==3,"CCZ injection inventory")
        require(all(ends[supplies[0]["id"]]<=starts[o["id"]] for o in injection),"CCZ input deadline")
    require(all(state[q]=="clean" for q in state if roles[q]=="injection_port"),"dirty final port")
    require(all(state[q]=="clean" or (state[q]=="measured" and spec.get("resource_reuse")=="accepted_replacement")
                for q in state if roles[q]=="CCZ_resource"),"unconsumed/uncleared resource")
    latency=max(ends,default=0)
    area=layout["rows"]*layout["cols"]
    require(result["latency_logical_cycles"]==latency,"latency total")
    require(result["reserved_patches"]==area and result["reserved_patch_cycles"]==area*latency,"area/volume total")
    require(result["primitive_counts"]==dict(counts) and result["operations"]==len(ops),"instruction inventory")
    require(result["final_patch_states"]==[state[q] for q in range(len(state))],"final lifecycle")
    return {"passed":True,"operations":len(ops),"latency_logical_cycles":latency,
            "reserved_patches":area,"reserved_patch_cycles":area*latency,
            "primitive_counts":dict(counts),"CCZ_states":counts["ccz_supply"],
            "CCZ_resource_sites":sum(r=="CCZ_resource" for r in roles.values()),
            "injection_port_sites":sum(r=="injection_port" for r in roles.values()),
            "module_times":{str(k):v for k,v in module_times.items()},
            "checks":["exact_input_inventory","quantum_and_classical_dependencies",
                      "branch_predicates","feedback_readiness","full_H_footprints",
                      "directed_full_CNOT_routes","space_time_exclusion",
                      "CCZ_supply_Bell_reset_lifecycle","module_homogeneous_stages",
                      "zero_duration_event_order","area_latency_volume"]}
