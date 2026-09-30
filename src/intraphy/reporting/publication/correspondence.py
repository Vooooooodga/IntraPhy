"""Whole-locus genomic correspondence and qualified exon-span observations."""
from .drawing import box, line, text, track, document, RULE, TEAL, GOLD, PURPLE, RED, MUTED


def draw():
    b = [text(42, 46, "Genomic correspondence and exon-span observations", 27, bold=True),
         text(42, 76, "Shared coordinates depend on qualified whole-locus sequence evidence", 15, color=MUTED)]
    # A: genomic DNA alignment is the primary correspondence evidence.
    b += [text(55, 131, "A   Genomic alignment", 18, bold=True),
          text(55, 160, "annotated physical exon intervals", 14, color=MUTED),
          text(55, 190, "Species 1", 13, color=MUTED)]
    b += track(62, 217, 500, 100, [(4,24),(31,45),(69,80),(87,98)], height=24)
    b += [text(55, 393, "Species 2", 13, color=MUTED)]
    b += track(62, 355, 500, 100, [(4,24),(31,45),(69,80),(87,98)], color=PURPLE, height=18)
    b += [line(182, 241, 182, 355, GOLD, 2), line(302,241,302,355,GOLD,2),
          line(412,241,412,355,GOLD,2), line(507,241,507,355,GOLD,2),
          text(61, 431, "Shared alignment coordinates", 13, color=MUTED),
          line(580, 280, 632, 280, arrow=True)]
    # B: interval matches are chained in local genomic order; no algorithmic
    # path-selection procedure is implied by this schematic.
    b += [text(650, 131, "B   Interval correspondence", 18, bold=True),
          text(650, 159, "collinear matches, ordered locally", 13, color=MUTED)]
    nodes = ((660,232,"Match 1"),(812,232,"Match 2"),(964,232,"Match 3"))
    for x,y,label in nodes:
        b += [box(x,y,132,46,"white",TEAL,5), text(x+66,y+31,label,13,anchor="middle")]
    b += [line(792,255,812,255,arrow=True),line(944,255,964,255,arrow=True),
          text(660, 324, "Order + strand + shared sequence", 13, color=MUTED),
          text(660, 374, "Sequence-qualified local chain", 13, color=MUTED),
          line(1138, 280, 1180, 280, arrow=True)]
    # C: model observations distinguish qualified exon spans, supported DNA
    # absence, and unknown/conflicting evidence.
    b += [text(1190, 131, "C   Tip observations", 19, bold=True)]
    rows = ((171, ("Qualified exon-span observation", "physical interval on genomic DNA"), TEAL, "E"),
            (268, ("DNA absence", "paired flanks support the call"), GOLD, "∅"),
            (365, ("Unknown or conflicting", "mapping / annotation unresolved"), RED, "?"))
    for y,labels,color,mark in rows:
        b += [box(1190,y,542,88,"white",RULE,5), text(1210,y+56,mark,23,color=color,bold=True),
              text(1250,y+39,labels[0],13), text(1250,y+75,labels[1],12)]
    b += [line(42, 473, 1732, 473, RULE, 1),
          text(48, 514, "Physical intervals define genomic exon-span geometry.", 14),
          text(48, 544, "Unknown status retains any known DNA constraints.", 14, color=MUTED)]
    return document(1780, 590, "Genomic correspondence and exon-span observations", b,
                    "Panel A emphasizes whole-locus genomic DNA alignment and shared coordinates. Panel B schematically shows collinear interval matches in local genomic order without specifying a path-selection algorithm. Panel C distinguishes a qualified physical exon-span observation, DNA absence supported by paired flanks, and unknown or conflicting evidence. Correspondence and observations remain conditional on sequence support and annotation coverage.")
