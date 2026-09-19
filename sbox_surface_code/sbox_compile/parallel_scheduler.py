"""Event list scheduling; no synthesis, placement changes, or new routes."""
import heapq
import random
import time


def prepare(model):
    nodes = model['nodes']
    index = {n['id']: i for i,n in enumerate(nodes)}
    successors = [[] for _ in nodes]
    indegree = [0]*len(nodes)
    for a,b in model['dependencies']:
        successors[index[a]].append(index[b]); indegree[index[b]] += 1
    coordinates = sorted({tuple(p) for n in nodes for p in n['footprint']})
    coordinate_index = {p:i for i,p in enumerate(coordinates)}
    masks = [sum(1 << coordinate_index[tuple(p)] for p in n['footprint']) for n in nodes]
    tails = [0]*len(nodes)
    for i in reversed(range(len(nodes))):
        tails[i] = nodes[i]['duration'] + max((tails[j] for j in successors[i]), default=0)
    group_of, acquisitions, releases, hold_masks = {}, {}, {}, []
    for k,g in enumerate(model['port_lifetimes']):
        a,b = index[g['acquire']],index[g['bell']]
        acquisitions[a]=k; releases[b]=k; group_of[a]=group_of[b]=k
        hold_masks.append(sum(1<<coordinate_index[tuple(p)] for p in g['ports']+g['resources']))
    return nodes,successors,indegree,masks,tails,group_of,acquisitions,releases,hold_masks


def run(prepared, seed=0, mode='critical', hint=None):
    nodes,succ,counts,masks,tails,group_of,acquire,release,hmasks = prepared
    pending=list(counts)
    ready={i for i,n in enumerate(pending) if not n}
    randomizer=random.Random(seed)
    noise=[randomizer.random() for _ in nodes]
    if mode=='hint':
        priority=lambda i:(hint[nodes[i]['id']]+(seed%16)*noise[i],-tails[i],i)
    elif mode=='canonical':
        priority=lambda i:(i,)
    elif mode=='small':
        priority=lambda i:(-tails[i],masks[i].bit_count(),i)
    elif mode=='jitter':
        priority=lambda i:(-(tails[i]+8*noise[i]),i)
    elif mode=='random':
        priority=lambda i:(-tails[i]*(.8+.4*noise[i]),i)
    else:
        priority=lambda i:(-tails[i],i)
    now=0; active=[]; held={}; events=[]; done=0; busy=0

    def complete(i):
        nonlocal done
        done+=1
        if i in release:
            del held[release[i]]
        for j in succ[i]:
            pending[j]-=1
            if pending[j]==0: ready.add(j)

    while done < len(nodes):
        while active and active[0][0]<=now:
            _,i=heapq.heappop(active)
            busy &= ~masks[i]
            complete(i)
        changed=True
        while changed:
            changed=False
            for i in sorted(ready,key=priority):
                forbidden=busy
                own=group_of.get(i)
                for k,mask in held.items():
                    if k!=own: forbidden |= mask
                if masks[i]&forbidden: continue
                duration=nodes[i]['duration']
                if i in acquire:
                    mask=hmasks[acquire[i]]
                    if mask&forbidden: continue
                    held[acquire[i]]=mask
                ready.remove(i)
                events.append(dict(id=nodes[i]['id'],start=now,end=now+duration))
                if duration:
                    busy |= masks[i]; heapq.heappush(active,(now+duration,i))
                else:
                    complete(i)
                changed=True
        if done==len(nodes): break
        if not active:
            raise RuntimeError('deadlock: unresolved dependency or occupied CCZ lifetime')
        now=active[0][0]
    return dict(events=events,latency=max(e['end'] for e in events),seed=seed,policy=mode,
                dependency_only_lower_bound=max(tails))


def compact(model,timeline):
    """Left-shift a feasible trace without reversing resource-use orders."""
    nodes=model['nodes']; index={n['id']:i for i,n in enumerate(nodes)}
    events={e['id']:e for e in timeline['events']}
    dispatch={e['id']:i for i,e in enumerate(timeline['events'])}
    edges=set(map(tuple,model['dependencies']))
    use={}
    own={}
    for g in model['port_lifetimes']:
        points=set(map(tuple,g['ports']+g['resources']))
        own[g['acquire']]=points; own[g['bell']]=points
        for p in points:
            use.setdefault(p,[]).append((events[g['acquire']]['start'],dispatch[g['acquire']],
                                         g['acquire'],g['bell']))
    for n in nodes:
        for p in set(map(tuple,n['footprint']))-own.get(n['id'],set()):
            use.setdefault(p,[]).append((events[n['id']]['start'],dispatch[n['id']],n['id'],n['id']))
    for entries in use.values():
        entries.sort()
        for left,right in zip(entries,entries[1:]):
            assert events[left[3]]['end']<=events[right[2]]['start'], 'invalid seed footprint order'
            if left[3]!=right[2]: edges.add((left[3],right[2]))
    succ=[[] for _ in nodes]; pending=[0]*len(nodes); earliest=[0]*len(nodes)
    for a,b in edges:
        succ[index[a]].append(index[b]); pending[index[b]]+=1
    ready=[i for i,n in enumerate(pending) if not n]; heapq.heapify(ready)
    result=[]; topological={}
    while ready:
        i=heapq.heappop(ready); n=nodes[i]; end=earliest[i]+n['duration']
        topological[n['id']]=len(result)
        result.append(dict(id=n['id'],start=earliest[i],end=end))
        for j in succ[i]:
            earliest[j]=max(earliest[j],end); pending[j]-=1
            if not pending[j]:heapq.heappush(ready,j)
    assert len(result)==len(nodes),'cyclic seed resource order'
    result.sort(key=lambda e:(e['start'],topological[e['id']]))
    return dict(events=result,latency=max(e['end'] for e in result),policy='resource_order_compaction')


def search(model, starts=64, initial=None):
    prepared=prepare(model)
    results=[]; best=None; cpu=time.process_time(); wall=time.monotonic()
    if initial is not None:
        best=compact(model,initial)
    hint={e['id']:e['start'] for e in best['events']} if best is not None else None
    for seed in range(starts):
        mode=['critical','canonical','small'][seed] if seed<3 else ('jitter' if seed%2 else 'random')
        if hint is not None and seed>=3 and seed%3:mode='hint'
        before=time.process_time()
        candidate=run(prepared,seed,mode,hint)
        raw_latency=candidate['latency']
        candidate=compact(model,candidate)
        results.append(dict(seed=seed,policy=mode,latency=candidate['latency'],raw_latency=raw_latency,
                            cpu_seconds=time.process_time()-before))
        if best is None or candidate['latency']<best['latency']:
            best=candidate
            best.update(seed=seed,initial_policy=mode)
            hint={e['id']:e['start'] for e in best['events']}
    return dict(best=best,runs=results,cpu_seconds=time.process_time()-cpu,
                elapsed_seconds=time.monotonic()-wall)
