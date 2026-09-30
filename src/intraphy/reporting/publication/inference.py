"""Default genomic exon-span CTMC pruning, fitting, and branch summaries."""
from .drawing import box, line, text, math_text, tree, track, document, PALE, RULE, TEAL, PURPLE, MUTED


def draw():
    b = [text(42, 50, "Default genomic exon-span CTMC on a supplied tree", 28, bold=True),
         text(42, 87, "One physical exon-span geometry with declared DNA-material states", 15, color=MUTED)]
    # A: local tip observations and supplied tree.
    b += [box(42,120,515,525,"white",RULE),
          text(66,163,"A   Tree and tip evidence",19,bold=True)]
    b += tree(86,235,.86)
    for i,s in enumerate("ABCD"):
        yy=235+i*60
        b.append(line(86+145*.86+42,yy,300,yy,RULE,1.5))
        b += track(300,yy,215,100,[(4,24),(37,59),(76,96)],height=18)
    b += [text(66, 510, "Physical exon-span compatibility", 14, color=MUTED),
          text(66, 548, "Root, topology and branch lengths", 14),
          text(66, 582, "are supplied model conditions", 14)]
    # B: pruning is conditional on one origin scenario; scenario likelihoods
    # are prior weighted before the linked-unit family composite likelihood.
    b += [box(590,120,1158,525,PALE,RULE),
          text(614,163,"B   Conditional pruning and family likelihood",19,bold=True)]
    tx,ty,scale=655,251,1.16
    b += tree(tx,ty,scale)
    xa,xb=tx+70*scale,tx+145*scale
    for yy,dy in ((ty+1,-7),(ty+60,8),(ty+120,-7),(ty+180,8)):
        b.append(line(xb-5,yy+dy,xa+8,yy+dy,PURPLE,2.5,arrow=True))
    for yy,dy in ((ty+30,-9),(ty+150,9)):
        b.append(line(xa-4,yy+dy,tx+8,yy+dy,PURPLE,2.5,arrow=True))
    b += [text(932,222,"Prune within one origin scenario s",14),
          math_text(932,267,[("L",False),("u",True),("(s)",False),(" = Σ",False),("i",True),(" π",False),("i",True),("(s)",False),(" L",False),("root",True),("(i | s)",False)],15,color=PURPLE,bold=True),
          text(932,315,"Scenario likelihoods are weighted by their prior",14),
          math_text(932,359,[("L",False),("u",True),(" = Σ",False),("s",True),(" p(s) L",False),("u",True),("(s)",False)],15,color=PURPLE,bold=True),
          text(932,407,"Family composite log likelihood",14,bold=True),
          math_text(932,450,[("ℓ",False),("family",True),(" = Σ",False),("u",True),(" log L",False),("u",True)],15,color=PURPLE,bold=True),
          math_text(932,504,[("P",False),("s,b",True),("(t) = exp(Q",False),("s,b",True),("t)",False)],14,color=PURPLE),
          text(932,544,"Node marginals integrate origin scenarios",14)]
    # C: default branch table summaries are endpoint contrasts.
    b += [box(42,690,530,260,"white",RULE),text(66,733,"C   Branch endpoint summaries",18,bold=True),
          text(72,784,"Parent",14,bold=True),text(385,784,"Child",14,bold=True)]
    b += track(72,810,185,100,[(5,38),(55,85)],height=20)
    b += [line(265,822,340,822,arrow=True),line(265,839,340,839,arrow=True)]
    b += track(353,810,185,100,[(5,38),(55,68),(77,85)],height=20)
    b += [text(66,884,"Branch P: geometry differs",12,color=PURPLE,bold=True),
          text(66,920,"Joint modal endpoint pair",12,color=MUTED)]
    # D: one scalar rate per family, shared across edit kinds and linked units.
    b += [box(630,690,530,260,PALE,RULE),text(654,733,"D   Per-family rate fit",18,bold=True),
          text(654,786,"One nonnegative μ per family",15,bold=True),
          text(654,824,"Shared across edit kinds",14),
          text(654,860,"and linked local units",14),
          text(654,903,"Linked units: composite likelihood",13,color=MUTED)]
    # E: conditions are explicit inputs/priors, not an automatic sensitivity run.
    b += [box(1218,690,530,260,"white",RULE),text(1242,733,"E   Declared conditions",18,bold=True),
          text(1242,786,"Supplied tree + branch lengths",14),
          text(1242,824,"Finite exon-boundary catalogue",14),
          text(1242,860,"Root and origin assumptions",14),
          text(1242,903,"Inference conditional on",12,color=MUTED),
          text(1242,933,"these declared choices",12,color=MUTED)]
    b += [math_text(48,990,[("q",False),("ij",True),(" = μ Σ",False),("e:i→j",True),(" w",False),("e",True),("   (i ≠ j)",False)],15,color=PURPLE,bold=True),
          math_text(48,1027,[("q",False),("ii",True),(" = −Σ",False),("j≠i",True),(" q",False),("ij",True)],15,color=PURPLE,bold=True),
          text(48,1063,"μ is shared across the edit graph; weights divide opportunities.",13),
          text(48,1098,"Equal-rate graph baseline ≠ biological rate equality; endpoints do not specify paths.",12,color=MUTED)]
    return document(1790, 1125, "Default genomic exon-span CTMC on a supplied tree", b,
                    "Panel A shows physical exon-span tip observations on a supplied rooted tree. Panel B shows pruning conditional on one declared material-origin scenario s, using the scenario-specific root prior π_i(s) and branch/origin-specific generator Q_s,b. Scenario likelihoods are weighted by their prior; linked-unit likelihoods contribute to the family composite log likelihood. Panel C shows branch endpoint geometry-change probability and joint modal parent/child configurations, not transition counts. Panel D shows one nonnegative rate per family shared across edit kinds and linked units. Panel E lists supplied tree and branch lengths, finite boundary catalogue, and root/origin assumptions. For eligible edit edges, q_ij = μ Σ_e w_e and q_ii = −Σ_j≠i q_ij. Equal rates define a graph baseline, not equal biological rates.")
