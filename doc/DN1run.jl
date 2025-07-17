using DN1, GraphRecipes, Graphs

"""
G = krožna_lestev(n)
Ustvari graf krožna lestev z `2n` točkami.
"""
function krožna_lestev(n)
    G = SimpleGraph(2 * n)
    # prvi cikel
    for i = 1:n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, 1, n)
    # drugi cikel
    for i = n+1:2n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, n + 1, 2n)
    # povezave med obema cikloma
    for i = 1:n
        add_edge!(G, i, i + n)
    end
    return G
end


"""
A = matrika(G::AbstractGraph, sprem)
Poišči matriko sistema linearnih enačb za vložitev grafa `G` s fizikalno metodo.
Argument `sprem` je vektor vozlišč grafa, ki nimajo določenih koordinat.
Indeksi v matriki `A` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function matrika(G, sprem)
    # preslikava med vozlišči in indeksi v matriki
    v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
    m = length(sprem)
    A = zeros(m, m)
    for i = 1:m
        vertex = sprem[i]
        sosedi = neighbors(G, vertex)
        for vertex2 in sosedi
            if haskey(v_to_i, vertex2)
                j = v_to_i[vertex2]
                A[i, j] = 1
            end
        end
        A[i, i] = -length(sosedi)
    end
    return A
end

"""
b = desne_strani(G::AbstractGraph, sprem, koordinate)
Poišči desne strani sistema linearnih enačb za eno koordinato vložitve grafa `G`
s fizikalno metodo. Argument `sprem` je vektor vozlišč grafa, ki nimajo
določenih koordinat. Argument `koordinate` vsebuje eno koordinato za vsa
vozlišča grafa. Metoda uporabi le koordinato vozlišč, ki so pritrjena.
Indeksi v vektorju `b` ustrezajo vozliščem v istem vrstnem redu,
kot nastopajo v argumentu `sprem`.
"""
function desne_strani(G, sprem, koordinate)
    set = Set(sprem)
    m = length(sprem)
    b = zeros(m)
    for i = 1:m
        v = sprem[i]
        for v2 in neighbors(G, v)
            if !(v2 in set) # dodamo le točke, ki so fiksirane
                b[i] -= koordinate[v2]
            end
        end
    end
    return b
end

"""
vloži!(G::AbstractGraph, fix, točke)
Poišči vložitev grafa `G` v prostor s fizikalno metodo. Argument `fix` vsebuje
vektor vozlišč grafa, ki imajo določene koordinate. Argument `točke` je
začetna vložitev grafa. Koordinate vozlišč, ki niso pritrjena, bodo nadomeščene
z novimi koordinatami.
Metoda ne vrne ničesar, ampak zapiše izračunane koordinate v matriko `točke`.
"""
function vloži!(G, fix, točke)
    sprem = setdiff(vertices(G), fix)
    dim, _ = size(točke)
    A = build_indexed_matrix(matrika(G, sprem))
    for k = 1:dim
        b = desne_strani(G, sprem, točke[k, :])
        x = conj_grad(A, b) 
        točke[k, sprem] = x
    end
end

m, n = 6, 6
G = Graphs.grid((m, n), periodic=false)
# vogali imajo stopnjo 2
vogali = filter(v -> degree(G, v) <= 2, vertices(G))
točke = zeros(2, n * m)
točke[:, vogali] = [0 0 1 1; 0 1 0 1]
vloži!(G, vogali, točke)
graphplot(G, x=točke[1, :], y=točke[2, :], curves=false)