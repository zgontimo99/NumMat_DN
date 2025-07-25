# Izračun ploščine znotraj zanke

Avtor: Timotej Zgonik

Paket, ki implementira osnovne računske operacije s polinomi in izračun integrala po Greenovem teoremu, demo skripta pa izvede izračun Bézierjeve krivulje, dane s kontrolnim poligonom: (0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (-1, 3), (0, 1), (1, 0). 

Skripto lahko poženemo z ukazom:
```jl
include("DN2/doc/DN2run.jl")
```
v interaktivni zanki Julie, potem ko smo ustrezno postavili paket v DN2.
## Testi
Teste poženemo z ukazom:
```
julia --project=DN2 -e "import Pkg; Pkg.test()"
```
## Poročilo PDF

Poročilo je bilo napisano s pomočjo spletnega urejevalnika Overleaf.
