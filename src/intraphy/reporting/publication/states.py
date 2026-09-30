"""Genomic exon-span geometries, material states, and tip compatibility."""
from .drawing import box, line, text, track, document, PALE, RULE, TEAL, PURPLE, GOLD, MUTED


def draw():
    b = [text(42, 46, "Genomic exon-span states and observation constraints", 27, bold=True),
         text(42, 77, "One ordered physical exon geometry per state", 15, color=MUTED)]
    # A: One observation unit has a single ordered physical exon-span geometry.
    b += [box(45, 120, 500, 280, PALE), text(68, 158, "A   Physical exon geometry", 19, bold=True),
          text(75, 217, "G₁", 16, bold=True)]
    b += track(126, 211, 340, 100, [(4,31),(49,72),(78,96)], height=27)
    b += [text(75, 282, "G₂", 16, bold=True)]
    b += track(126, 276, 340, 100, [(4,31),(49,96)], color=PURPLE, height=27)
    b += [text(68, 338, "Ordered physical intervals", 14),
          text(68, 375, "including terminal and UTR spans", 12, color=MUTED)]
    # B: Declared material tracts have a single root-or-branch origin.
    b += [box(582, 120, 520, 280, "white", RULE), text(605, 158, "B   Declared DNA-tract states", 19, bold=True),
          box(606, 214, 112, 48, PALE, RULE), text(662, 245, "0", 17, anchor="middle", bold=True),
          box(786, 214, 112, 48, PALE, TEAL), text(842, 245, "1", 17, anchor="middle", bold=True),
          box(966, 214, 112, 48, PALE, PURPLE), text(1022, 245, "2", 17, anchor="middle", bold=True),
          line(718,238,786,238,arrow=True), line(898,238,966,238,arrow=True),
          text(662, 294, "unintroduced", 10, anchor="middle", color=MUTED),
          text(842, 294, "present", 10, anchor="middle", color=TEAL),
          text(1022, 294, "deleted", 10, anchor="middle", color=PURPLE),
          text(608, 330, "One origin: root or one branch", 12, bold=True),
          text(608, 375, "Absent: 0 or 2; present: 1", 12, color=MUTED)]
    # C: Unknown exon status leaves geometry unconstrained subject to known DNA.
    b += [box(1140, 120, 602, 280, PALE), text(1163, 158, "C   Tip compatibility", 19, bold=True),
          text(1163, 198, "Unknown exon call; DNA present", 12, color=MUTED)]
    for x,label,detail,weight in ((1210,"G₁","m=1","w=1"),(1360,"G₂","m=1","w=1"),(1510,"G₃","m=0/2","w=0")):
        color = TEAL if weight == "w=1" else GOLD
        b += [f'<circle cx="{x}" cy="236" r="24" fill="white" stroke="{color}" stroke-width="2"/>',
              text(x,242,label,12,anchor="middle",bold=True),
              text(x,294,detail,10,anchor="middle",color=color),
              text(x,326,weight,10,anchor="middle",color=color,bold=True)]
    b += [text(1163, 361, "Compatible geometries get weight 1", 11),
          text(1163, 394, "Known DNA constraints remain active", 11, color=MUTED)]
    # D: Candidate boundaries and material cuts define the complete finite set.
    b += [line(42, 411, 1742, 411, RULE, 1), text(48, 453, "D   Finite geometry catalogue", 18, bold=True),
          text(48, 482, "Observed boundaries + tract cuts", 12, color=MUTED),
          line(545, 465, 640, 465, arrow=True), text(670,453,"Finite catalogue",14,bold=True),
          text(670,482,"Complete enumeration",12,color=MUTED)]
    for x,y,label in ((1015,450,"G₁"),(1155,450,"G₂"),(1085,520,"G₃"),(1235,520,"…")):
        b += [f'<circle cx="{x}" cy="{y}" r="20" fill="{PALE}" stroke="{PURPLE}" stroke-width="2"/>',
              text(x,y+6,label,12,anchor="middle",bold=True)]
    b += [line(1035,450,1135,450,arrow=True),line(1168,466,1100,504,arrow=True),
          line(1105,520,1215,520,arrow=True),
          text(1320, 468, "No default count cap", 14, bold=True),
          text(1320, 506, "Limits stop; no truncation", 11)]
    return document(1790, 570, "Genomic exon-span states and observation constraints", b,
                    "Panel A shows example single physical exon-span geometries as separate finite states, including terminal and UTR spans. Panel B shows declared DNA material-tract states 0 (unintroduced), 1 (present), and 2 (deleted), with one origin at the root or on one branch; tip absence is compatible with 0 or 2, while presence requires 1. Panel C shows unknown exon status with compatibility weight one for geometries satisfying recorded DNA constraints. Panel D shows the complete finite state space conditional on observed exon boundaries and declared tract cuts; there is no default state-count cap.")
