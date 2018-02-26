(TeX-add-style-hook
 "Research_Statement"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "Geoff"
    "graphicx"
    "subcaption")
   (LaTeX-add-labels
    "sec:intro"
    "sec:method"
    "sec:experiment"
    "sec:futurework")
   (LaTeX-add-bibliographies
    "../../BibTex/research"))
 :latex)

