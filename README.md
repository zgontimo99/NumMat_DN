# Izračun periode limitnega cikla Van der Polovega oscilatorja

Avtor: Timotej Zgonik

Paket, ki implementira metodo Runge-Kutta četrtega reda in definicijo Van der Polovega oscilatorja, demo skripta pa izvede izračun periode limitnega cikla oscilatorja pri μ = 4.
Skripto lahko poženemo z ukazom:
```jl
include("DN3/doc/DN3run.jl")
```
v interaktivni zanki Julie, potem ko smo ustrezno postavili paket v DN3.
## Testi
Teste poženemo z ukazom:
```
julia --project=DN3 -e "import Pkg; Pkg.test()"
```
## Poročilo PDF

Poročilo je bilo napisano s pomočjo spletnega urejevalnika Overleaf.
