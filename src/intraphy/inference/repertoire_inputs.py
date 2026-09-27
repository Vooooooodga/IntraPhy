"""Strict versioned JSON I/O for conditional repertoire CTMC models."""
from __future__ import annotations
from dataclasses import dataclass, field
import json, math
from pathlib import Path
import numpy as np

from ..structure.types import Catalogue, ExonConfiguration, ExonSpan
from ..structure.space import StateSpace
from ..structure.repertoire import ExonRepertoire, repertoire_compatibility
from ..structure.repertoire_process import RepertoireProcess, RepertoireEvent, deletion_process
from ..structure.repertoire_edits import coupled_process
from .repertoire_rates import RepertoireInferenceUnit

SCHEMA = "intraphy.exon-repertoire-model/1"

@dataclass(frozen=True)
class RepertoireBundle:
    units: tuple[RepertoireInferenceUnit, ...]
    rates: dict[tuple[str, str], dict[str, float]]
    rate_provenance: dict[tuple[str, str], str]
    branch_length_unit: str
    scale: float
    provenance: str
    observation_metadata: dict[tuple[str, str], dict]
    model_record: dict = field(default_factory=dict)

def _obj(x, keys, name):
    if not isinstance(x, dict) or set(x) != set(keys): raise ValueError(f"strict {name} fields required")
    return x
def _span(x):
    if not isinstance(x, list) or len(x) != 2 or any(type(v) is not int for v in x): raise ValueError("invalid exon span")
    return ExonSpan(*x)
def _rep(x):
    _obj(x, {"material", "configurations"}, "repertoire")
    if not isinstance(x["material"], list) or any(type(v) is not int or v not in (0,1,2) for v in x["material"]): raise ValueError("invalid material")
    if not isinstance(x["configurations"], list): raise ValueError("configurations must be a list")
    return ExonRepertoire(tuple(_config(c) for c in x["configurations"]), tuple(x["material"]))
def _config(x):
    _obj(x,{"exons","material"},"configuration")
    if not isinstance(x["exons"], list) or not isinstance(x["material"], list) or any(type(v) is not int or v not in (0,1,2) for v in x["material"]): raise ValueError("invalid configuration containers")
    return ExonConfiguration(tuple(_span(e) for e in x["exons"]),tuple(x["material"]))
def encode_repertoire(r):
    return {"material": list(r.material), "configurations":[{"exons":[[e.start,e.end] for e in c.exons],"material":list(c.material)} for c in r.configurations]}
def _finite(v, name):
    if isinstance(v, bool) or not isinstance(v,(int,float)): raise ValueError(f"invalid {name}")
    try: value=float(v)
    except (OverflowError, ValueError): raise ValueError(f"invalid {name}")
    if not math.isfinite(value): raise ValueError(f"invalid {name}")
    return value
def _strict_load(path):
    def pairs(pairs):
        d={}
        for k,v in pairs:
            if k in d: raise ValueError(f"duplicate JSON key: {k}")
            d[k]=v
        return d
    def finite_float(x):
        v=float(x)
        if not math.isfinite(v): raise ValueError("nonfinite JSON number")
        return v
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON number")), parse_float=finite_float)
def read_repertoire_model(path, catalogues, tree):
    data=_strict_load(path); _obj(data,{"schema","model","branch_length_unit","scale","provenance","units"},"model")
    if data["schema"]!=SCHEMA or data["model"]!="exon-repertoire-ctmc": raise ValueError("schema/model mismatch")
    branch=data["branch_length_unit"]; provenance=data["provenance"]
    if not isinstance(branch,str) or not branch.strip() or not isinstance(provenance,str) or not provenance.strip(): raise ValueError("provenance and branch unit required")
    scale=_finite(data["scale"],"scale")
    if scale<=0 or not isinstance(data["units"],list) or not data["units"]: raise ValueError("invalid units")
    from ..structure.validation import validate_collection
    catalogues=tuple(validate_collection(tuple(catalogues)))
    cats={(c.family,c.unit):c for c in catalogues}; units=[]; rates={}; rprov={}; ometa={}; seen=set()
    leaves=set(tree.leaf_by_label)
    for _, child in tree.edges():
        length=tree.branch_length(child)
        if not isinstance(length,(int,float)) or isinstance(length,bool) or not math.isfinite(length) or length<0: raise ValueError("invalid tree branch length")
    for u in data["units"]:
        _obj(u,{"family","unit","process","rates","rate_provenance","root","observation"},"unit")
        if not isinstance(u["family"],str) or not u["family"].strip() or not isinstance(u["unit"],str) or not u["unit"].strip(): raise ValueError("invalid family/unit")
        key=(u["family"],u["unit"])
        if key in seen or key not in cats: raise ValueError("catalogue unit mismatch or duplicate")
        seen.add(key); c=cats[key]; process_data=u["process"]; pprov=process_data.get("provenance") if isinstance(process_data,dict) else None
        if any(o.species not in leaves for o in c.observations): raise ValueError("catalogue observation species is outside tree")
        mode=process_data.get("mode") if isinstance(process_data,dict) else None
        if not isinstance(mode,str): raise ValueError("process mode must be a string")
        if mode not in {"explicit","deletion","coupled"}: raise ValueError("unsupported process mode")
        if mode=="explicit":
            _obj(process_data,{"mode","states","events","provenance"},"explicit process")
            if not isinstance(process_data["states"],list) or not isinstance(process_data["events"],list): raise ValueError("process states/events must be lists")
            states=tuple(_rep(x) for x in process_data["states"]); events=[]
            for e in process_data["events"]:
                _obj(e,{"source","target","kind","opportunity","weight"},"event")
                if type(e["source"]) is not int or type(e["target"]) is not int or not 0<=e["source"]<len(states) or not 0<=e["target"]<len(states): raise ValueError("event index")
                events.append(RepertoireEvent(states[e["source"]],states[e["target"]],e["kind"],e["opportunity"],_finite(e["weight"],"weight")))
            process=RepertoireProcess(states,tuple(events),pprov,c)
        elif mode=="deletion":
            if set(process_data) - {"mode","seeds","provenance","max_states"} or not {"mode","seeds","provenance"} <= set(process_data): raise ValueError("invalid deletion process fields")
            if not isinstance(process_data["seeds"],list): raise ValueError("seeds must be a list")
            process=deletion_process(c,tuple(_rep(x) for x in process_data["seeds"]),provenance=pprov,max_states=process_data.get("max_states"))
        else:
            required={"mode","policy","option_bank","insertion_outcomes","seeds","provenance"}
            if set(process_data)-required-{"max_states"} or not required <= set(process_data): raise ValueError("invalid coupled process fields")
            def event(x):
                _obj(x,{"source","target","kind","opportunity","weight"},"insertion event")
                return RepertoireEvent(_rep(x["source"]),_rep(x["target"]),x["kind"],x["opportunity"],_finite(x["weight"],"weight"))
            if not all(isinstance(process_data[k],list) for k in ("seeds","option_bank","insertion_outcomes")): raise ValueError("coupled process collections must be lists")
            process=coupled_process(c,tuple(_rep(x) for x in process_data["seeds"]),policy=process_data["policy"],option_bank=tuple(_config(x) for x in process_data["option_bank"]),insertion_outcomes=tuple(event(x) for x in process_data["insertion_outcomes"]),provenance=pprov,max_states=process_data.get("max_states"))
        root=u["root"]; _obj(root,{"entries","provenance"},"root")
        if not isinstance(root["provenance"],str) or not root["provenance"].strip(): raise ValueError("root provenance required")
        if not isinstance(root["entries"],list): raise ValueError("root entries must be a list")
        process_index=process.index
        prior=np.zeros(len(process.states)); used=set()
        for ent in root["entries"]:
            _obj(ent,{"state","weight"},"root entry"); r=_rep(ent["state"]); idx=process_index.get(r)
            if idx is None or idx in used: raise ValueError("root state mismatch or duplicate")
            used.add(idx); prior[idx]=_finite(ent["weight"],"root weight")
            if prior[idx] < 0: raise ValueError("root weights must be nonnegative")
        if not np.isclose(prior.sum(),1.,rtol=0,atol=1e-12): raise ValueError("root prior must sum to one")
        obs=u["observation"]
        if not isinstance(obs,dict): raise ValueError("observation must be an object")
        required={"mode","provenance","tips"} if obs.get("mode")=="explicit" else {"mode","provenance","view","missing_policy","complete_species"}
        if set(obs)!=required: raise ValueError("strict observation fields required")
        if not isinstance(obs["provenance"],str) or not obs["provenance"].strip(): raise ValueError("observation provenance required")
        if not isinstance(obs["mode"],str): raise ValueError("observation mode must be a string")
        if obs["mode"]=="explicit":
            if not isinstance(obs["tips"],dict): raise ValueError("tips must be an object")
            if set(obs["tips"])!=leaves: raise ValueError("explicit tips must cover every leaf")
            if any(not isinstance(v,list) for v in obs["tips"].values()): raise ValueError("invalid tip emission")
            tips={k:np.asarray([_finite(x,"tip emission") for x in v],dtype=float) for k,v in obs["tips"].items()}
            if any(v.shape!=(len(process.states),) or not np.isfinite(v).all() or (v<0).any() or (v>1).any() for v in tips.values()): raise ValueError("invalid tip emission")
        elif obs["mode"]=="constraints":
            if not isinstance(obs["view"],str) or not isinstance(obs["missing_policy"],str) or obs["view"] not in {"annotation","evidence"} or obs["missing_policy"] not in {"uninformative"}: raise ValueError("invalid constraint observation policy")
            if not isinstance(obs["complete_species"],list) or any(not isinstance(s,str) or not s.strip() for s in obs["complete_species"]) or len(set(obs["complete_species"])) != len(obs["complete_species"]): raise ValueError("invalid complete_species")
            if any(s not in leaves for s in obs["complete_species"]): raise ValueError("invalid complete_species")
            by_species={o.species:o for o in c.observations}
            if any(s not in by_species or by_species[s].kind in {"unknown","partial","excluded"} for s in obs["complete_species"]): raise ValueError("complete species must have resolved observations")
            members=tuple(sorted({m for s in process.states for m in s.configurations})); ss=StateSpace(c,members,(),False,"declared_repertoire_projection",len(members))
            tips={}
            for leaf in sorted(leaves):
                o=by_species.get(leaf)
                tips[leaf]=np.ones(len(process.states)) if o is None else repertoire_compatibility(ss,process.states,o,obs["view"],complete_repertoire=(leaf in obs["complete_species"]))
        else: raise ValueError("observation mode")
        if not isinstance(u["rates"],dict) or not isinstance(u["rate_provenance"],str): raise ValueError("invalid rate record")
        rates[key]={k:_finite(v,"rate") for k,v in u["rates"].items()}; rprov[key]=u["rate_provenance"]; ometa[key]=obs
        kinds={e.kind for e in process.events}
        if set(rates[key]) != kinds or any(v<0 for v in rates[key].values()) or not isinstance(rprov[key],str) or not rprov[key].strip(): raise ValueError("invalid rate record")
        units.append(RepertoireInferenceUnit(process,tree,tips,prior,root["provenance"],obs["provenance"]))
    if seen != set(cats): raise ValueError("catalogue collection mismatch")
    return RepertoireBundle(tuple(units),rates,rprov,branch,scale,provenance,ometa,data)
