"""A two-row overview of the default genomic exon-span CTMC."""
from .drawing import box, line, text, tree, track, document, RULE, TEAL, MUTED, GOLD, PURPLE, PALE


def panel(x, y, title, number):
    return [box(x, y, 540, 455, "white", RULE),
            text(x+24, y+43, number, 18, color=TEAL, bold=True),
            text(x+65, y+43, title, 19, bold=True)]


def draw():
    b = [text(42, 54, "Genomic exon-span evolution on a supplied tree", 28, bold=True),
         text(42, 88, "Physical exon geometry, declared DNA material, and conditional phylogenetic summaries", 15, color=MUTED)]
    cells = ((42,120),(625,120),(1208,120),(42,625),(625,625),(1208,625))
    titles = ("Orthologous loci", "Genomic correspondence", "Physical exon spans",
              "Admissible edits", "Likelihood inference", "Branch summaries")
    for i, ((x,y), title) in enumerate(zip(cells,titles),1):
        b += panel(x,y,title,"ABCDEF"[i-1])
    # A: Input loci and the supplied species tree.
    b += tree(80, 233, .75, labels=False)
    for i,s in enumerate("ABCD"):
        yy=233+i*60
        b.append(line(189,yy,241,yy,RULE,1.5))
        b += [text(247,yy+7,s,14,bold=True)]
        b += track(276,yy,240,100,[(4,22),(35,48),(66,80),(88,98)],height=17)
    b += [text(68, 511, "Genome + exon annotation", 14, color=MUTED),
          text(68, 544, "Rooted species tree supplied", 14, color=MUTED)]
    # B: Whole-locus genomic alignment defines shared comparison coordinates.
    x,y=cells[1]
    b += [text(x+30,y+105,"Whole-locus genomic alignment",14,color=MUTED)]
    b += track(x+30,y+133,420,100,[(3,24),(29,50),(67,83),(88,98)],height=24)
    b += [line(x+150,y+157,x+150,y+267,GOLD,2),line(x+294,y+157,x+294,y+267,GOLD,2),
          line(x+399,y+157,x+399,y+267,GOLD,2),text(x+30,y+294,"Shared alignment coordinates",14,color=MUTED)]
    b += track(x+30,y+319,420,100,[(3,24),(29,50),(67,83),(88,98)],color=PURPLE,height=16)
    b += [text(x+30,y+438,"Sequence · strand · local order",13,color=MUTED)]
    # C: One physical interval geometry per species; UTR and terminal spans remain.
    x,y=cells[2]
    for i,(s,exons) in enumerate((("A",[(3,26),(39,61),(74,97)]),
                                  ("B",[(3,26),(39,97)]),
                                  ("C",[(3,26),(56,61),(74,97)]),
                                  ("D",[(3,26),(74,97)]))):
        yy=y+153+i*61
        b += [text(x+30,yy+6,s,14,bold=True)] + track(x+70,yy,405,100,exons,height=19)
    b += [text(x+30,y+400,"Physical intervals incl. UTR spans",13,color=MUTED),
          text(x+30,y+434,"One geometry per observation unit",13,color=MUTED)]
    # D: The finite graph contains exon-span and declared-material edits.
    x,y=cells[3]
    for cx,cy,label in ((x+145,y+190,"C₁"),(x+365,y+190,"C₂"),(x+255,y+330,"C₃")):
        b += [f'<circle cx="{cx}" cy="{cy}" r="47" fill="{PALE}" stroke="{TEAL}" stroke-width="3"/>',
              text(cx,cy+8,label,18,anchor="middle",bold=True)]
    b += [line(x+192,y+190,x+318,y+190,arrow=True),line(x+354,y+230,x+288,y+288,arrow=True),
          text(x+35,y+405,"Exon and DNA edits",13,color=MUTED)]
    # E: The default inference path is one genomic exon-span CTMC.
    x,y=cells[4]
    b += tree(x+63,y+135,.78)
    b += [text(x+270,y+158,"Exon-span CTMC",13,bold=True),
          text(x+270,y+194,"ML fitting",13,color=PURPLE),
          text(x+270,y+230,"Tree pruning",13,color=PURPLE),
          text(x+270,y+266,"One family rate",13,bold=True),
          text(x+35,y+404,"Supplied rooted species tree",13,color=MUTED),
          text(x+35,y+438,"Shared rate across edit kinds",13,color=MUTED)]
    # F: Default branch output is endpoint-based, not a transition count.
    x,y=cells[5]
    b += tree(x+45,y+145,.66)
    b += [text(x+78,y+160,"pᵢ",13,anchor="middle",color=PURPLE,bold=True),
          text(x+35,y+371,"Ancestral-state marginals",13,bold=True),
          text(x+35,y+400,"Branch: P(geometry differs)",13,color=MUTED),
          text(x+35,y+429,"Joint modal parent → child",13,color=MUTED)]
    # Arrows show reading order within rows.
    for ax,ay,bx,by in ((582,350,612,350),(1165,350,1195,350),
                        (582,850,612,850),
                        (1165,850,1195,850)):
        b.append(line(ax,ay,bx,by,TEAL,2,arrow=True))
    return document(1790, 1110, "Genomic exon-span evolution on a supplied tree", b,
                    "Panels A–F summarize the default genomic exon-span CTMC: whole-locus genomic alignment, physical exon intervals including UTR spans, admissible exon/material edits, sum-product pruning with one fitted rate per family, ancestral-state marginals, and branch endpoint-change probabilities with joint modal configurations. Branch summaries describe endpoints, not transition counts or molecular mechanisms.")
