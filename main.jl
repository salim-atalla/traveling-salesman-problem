
#= 
    ATALLA - Salim
    ABDULLA - Rahaf
=#

using JuMP, GLPK



# ====================================================
# =========== Modèle Miller-Tucker-Zemlin ============
# ====================================================

#=
    Résolution exacte en utilisant le modèle Miller-Tucker-Zemlin
=#
function methodeMTZ(C::Matrix{Int64})

    # Nombre de lieux à visiter
    n::Int64 = size(C, 1)

    # Créer un modèle vide
    m::Model = Model(GLPK.Optimizer)

    # Déclarer les variables de décision
    @variable(m, x[1:n, 1:n], binary = true) # xij : 1 si le drône se rend directement du lieu i au lieu j, 0 sinon
    @variable(m, t[2:n]) # tj : Date à laquelle la ville j est visitée.

    # Déclaration de la fonction objectif
    @objective(m, Min, sum(sum(C[i,j]x[i,j] for j in 1:n) for i in 1:n))

    # # Déclaration des contraintes
    # Les contraintes des regroupements
    # @constraint(m, Contr1[i in 1:n], sum(x[i,j] for j in 1:n if i != j) == 1)
    for i in 1:n
        aff = @expression(m, 0.0x[i,i])
        for j in 1:n
            if i != j
                add_to_expression!(aff, 1.0, x[i,j])
            end
        end
        @constraint(m, aff == 1)
    end

    # @constraint(m, Contr2[j in 1:n], sum(x[i,j] for i in 1:n if i != j) == 1)
    for j in 1:n
        aff = @expression(m, 0.0x[j,j])
        for i in 1:n
            if i != j
                add_to_expression!(aff, 1.0, x[i,j])
            end
        end
        @constraint(m, aff == 1)
    end

    # La contrainte sur les dates
    # @constraint(m, Contr3[i=2:n,j=2:n], (t[i] - t[j] + n*x[i,j]) <= (n - 1))
    for i in 2:n, j in 2:n
        aff = @expression(m, n*x[i,j])
        add_to_expression!(aff, -1.0, t[i])
        add_to_expression!(aff, 1.0, t[j])
        @constraint(m, aff <= n - 1)
    end

    return m
end





# ====================================================
# ======== Modèle Dantzig-Fulkerson-Johnson ==========
# ====================================================

#=
    Récupérer le lieu le plus proche pour chaque lieu à partir d'une matrice xij
    donnée (d'après la résolution des PL)
=#
function permutation(x::Matrix{Float64})

    # Nombre de lignes/colonnes
    n::Int64 = size(x, 1)
    # vecteur ou on va stocker la résultat
    res::Vector{Int64} = Vector{Int64}(undef, n) 

    for i in 1:n
        j = 1
        while x[i,j] != 1
            j = j + 1
        end
        res[i] = j
    end
    return res
end




#=
    Créer une structure de donnée du type vecteur de vecteur d'entiers 'Vector{Vector{Int64}}' pour stocker ls cycles,
    les valeurs sont des vecteurs qui contiennent les lieux passés dans chaque cycle
=#
function creer_cycles(V::Vector{Int64})

    # La taille du vecteur
    n::Int64 = length(V)
    cycles::Vector{Vector{Int64}} = []
    ajoutee::Vector{Bool}  = falses(n) # true si le lieu est déja dans un cycle, false sinon

    for lieu in 1:n
        if !(ajoutee[lieu]) # vérifier que le lieu n'est pas déjà dans un cycle

            # Ajouter le premier lieu dans le cycle
            push!(cycles, [lieu]) 
            ajoutee[lieu] = true
            comp = V[lieu] 

            i = 1
            while comp != lieu && i <= n # chercher la fin de chaque cycle, et ajouter les lieux passés
                push!(cycles[end], comp)
                ajoutee[comp] = true
                comp = V[comp]
                i = i + 1
            end
        end
    end

    return cycles
end




#=
    Afficher les cycles
=#
function afficher_cycles(cycles::Vector{Vector{Int64}})

    for cycle in cycles
        print("(")
        for lieu in cycle
            print(lieu)
        end
        print(")")
    end
    println()
end





#=
    Retourner le cycle à casser (le plus petit cycle) à partir d'un dictionnaire des cycles
=#
function cycle_a_casser(cycles::Vector{Vector{Int64}})

    # chercher le cycle le plus petit
    min_len::Int64 = 1000000
    min_vec::Vector{Int64} = cycles[1]

    for cycle in cycles
        if length(cycle) < min_len
            min_len = length(cycle)
            min_vec = cycle
        end
    end

    return min_vec
end





#=
    Générer l'expression affine
=#
function generer_aff(m::Model, x::Array{VariableRef}, vec::Vector{Int64})

    n::Int64 = size(x, 1)
    aff::AffExpr = @expression(m, 0.0)
    for i in 1:n
        for j in 1:n
            if i != j && (i in vec) && (j in vec) 
                add_to_expression!(aff, 1.0, x[i,j])
            end
        end
    end
    return aff
end





#=
    Résolution exacte en utilisant le modèle Dantzig-Fulkerson-Johnson
=#
function methodeDFJ(C::Matrix{Int64})

    # Nombre de lieux à visiter
    n::Int64 = size(C, 1)

    # Créer un modèle vide
    m::Model = Model(GLPK.Optimizer)

    # Déclarer les variables de décision
    @variable(m, x[1:n, 1:n], binary = true) # xij : 1 si le drône se rend directement du lieu i au lieu j, 0 sinon

    # Déclaration de la fonction objectif
    @objective(m, Min, sum(sum(C[i,j]x[i,j] for j in 1:n) for i in 1:n))

    # Déclaration des contraintes
    # Les contraintes des regroupements
    # @constraint(m, Contr1[i in 1:n], sum(x[i,j] for j in 1:n if i != j) == 1)
    for i in 1:n
        aff = @expression(m, 0.0x[i,i])
        for j in 1:n
            if i != j
                add_to_expression!(aff, 1.0, x[i,j])
            end
        end
        @constraint(m, aff == 1)
    end

    # @constraint(m, Contr2[j in 1:n], sum(x[i,j] for i in 1:n if i != j) == 1)
    for j in 1:n
        aff = @expression(m, 0.0x[j,j])
        for i in 1:n
            if i != j
                add_to_expression!(aff, 1.0, x[i,j])
            end
        end
        @constraint(m, aff == 1)
    end
    
    # Résoudre le modèle pour la première fois
    optimize!(m)

    # --- Générer les contraintes de cassage des cycles ---
    cycles::Vector{Vector{Int64}} = creer_cycles(permutation(value.(x)))

    # nbCycles = 1
    while length(cycles) > 1

        # # Afficher les cycles
        # print("cycle ", nbCycles, ":\t")
        # afficher_cycles(cycles)
        # nbCycles = nbCycles + 1

        vec::Vector{Int64} = cycle_a_casser(cycles)
        aff::AffExpr = generer_aff(m, x, cycles[findfirst(c -> c[1] == vec[1], cycles)])
        @constraint(m, aff <= length(vec) - 1)

        # # Afficher le cycle à casser
        # println("CYCLE A CASSER: ", vec)

        optimize!(m)
        cycles = creer_cycles(permutation(value.(x)))
    end

    # # Afficher le cycle final
    # print("cycle ", nbCycles, ":\t")
    # afficher_cycles(cycles)

    return m
end





# ====================================================
# ================ Algorithme Glouton ================
# ====================================================

#=
    Résolution par un algorithme glouton (Algorithme itérative)
=#
function methodeGlouton(C::Matrix)

    # Nombre de lieux à visiter
    n::Int64 = size(C, 1)
    visite::Vector{Bool}  = falses(n)  # visite[i] est vrai si le lieu i a déjà été visité, false sinon
    solution::Vector{Int64} = Vector{Int64}(undef, n+1) # vecteur qui stocke l'ordre des lieux visités
    distance::Int64 = 0 # Distance totale entre les lieux

    # On commence au lieu 1, qui est l'origine et la destination
    solution[1] = 1 
    visite[1] = true
    
    for i in 2:n
        distance_min = 1000000 # grande valeur pour assurer distance_min est initialiser au premier lieu visité
        lieu_plus_proche = 0

        # On recherche le lieu non visité le plus proche du lieu courant
        for j in 1:n
            if !visite[j] && C[solution[i-1],j] < distance_min
                distance_min = C[solution[i-1],j]
                lieu_plus_proche = j
            end
        end
        solution[i] = lieu_plus_proche
        visite[lieu_plus_proche] = true
        distance += distance_min
    end
    
    # On retourne la solution, en ajoutant le lieu de départ à la fin pour compléter le cycle
    solution[end] = 1
    # Ajouter la distance entre le dernier lieu visité et le premier lieu visité à la distance totale
    distance += C[solution[n],1]
    
    return (solution, distance)
end





# ====================================================
# ================== Résolution TSP ==================
# ====================================================

#=
    Une fonction qu'on devra écrire (en fait, il y en aura une pour chaque méthode, plus d'autres fonctions utiles...)
=#
function resolutionRapideDuTSP(C::Matrix{Int64})

    # <--- Résolution avec les algorithmes de résolution exacte --->
    # # Modèle Miller-Tucker-Zemlin
    # m::Model = methodeMTZ(C)

    # # Modèle Dantzig-Fulkerson-Johnson
    m::Model = methodeDFJ(C)


    # Résolution (pour les méthodes MTZ et DFJ)
    optimize!(m)

    # Affichage des résultats
    status = termination_status(m)

    if status == MOI.OPTIMAL
        println("Problème résolu à l'optimalité")
        println("z = ", objective_value(m)) # affichage de la valeur optimale
        # println("x = ", value.(m[:x]))

    elseif status == MOI.INFEASIBLE
        println("Problème impossible!")

    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Problème non borné!")

    end


    # <--- Résolution avec l'algorithme de résolution aprochée --->

    # # Algorithme glouton
    # solution::Vector{Int64}, distance::Int64 = methodeGlouton(C)
    # println("La solution est : ", solution)
    # println("La distance totle entre les lieux est : ", distance)

end





#=
    fonction qui prend en paramètre un fichier contenant un distancier et qui retourne le tableau bidimensionnel correspondant
=#
function parseTSP(nomFichier::String)

    # Ouverture d'un fichier en lecture
    f::IOStream = open(nomFichier,"r")

    # Lecture de la première ligne pour connaître la taille n du problème
    s::String = readline(f) # lecture d'une ligne et stockage dans une chaîne de caractères
    tab::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false)) # Segmentation de la ligne en plusieurs entiers, à stocker dans un tableau (qui ne contient ici qu'un entier)
    n::Int64 = tab[1]

    # Allocation mémoire pour le distancier
    C = Matrix{Int64}(undef,n,n)

    # Lecture du distancier
    for i in 1:n
        s = readline(f)
        tab = parse.(Int64,split(s," ",keepempty = false))
        for j in 1:n
            C[i,j] = tab[j]
        end
    end

    # Fermeture du fichier
    close(f)

    # Retour de la matrice de coûts
    return C
end





#= 
    Exemple de script, qui résout ici les instances jusqu'à une taille de 40. Ce script devra être modifié/adapté suivant les besoins. 
    La macro @time mesure le temps (elapsed) d'exécution d'une fonction
    La consommation mémoire est aussi indiquée mais la valeur indiquée n'est correcte que pour un code 100% codé en Julia (et GLPK est codé en C) 
=#
function ExempleScriptTSP()

    # Première exécution sur l'exemple pour forcer la compilation si elle n'a pas encore été exécutée
    C::Matrix{Int64} = parseTSP("plat/exemple.dat")
    resolutionRapideDuTSP(C)

    file::String = ""

    println("========================================")
    println("================= Plat =================")
    println("========================================")

    # Série d'exécution avec mesure du temps pour des instances symétriques
    for i in 10:10:150
        file = "plat/plat$i.dat"
        C = parseTSP(file)
        println("Instance à résoudre : plat$i.dat")
        @time resolutionRapideDuTSP(C)
        println("-----------------------------")
    end

    println()
    println("========================================")
    println("================ Relief ================")
    println("========================================")

    # Série d'exécution avec mesure du temps pour des instances asymétriques
    for i in 10:10:150
        file = "relief/relief$i.dat"
        println("Instance à résoudre : relief$i.dat")
        C = parseTSP(file)
        @time resolutionRapideDuTSP(C)
        println("-----------------------------")
    end
end

# # Résoudre tous les instances
# ExempleScriptTSP()




#=
    Pour tester une seule instance
=#
function resoudre_une_instance(chemin::String) 
    @time resolutionRapideDuTSP(parseTSP(chemin))
end

# resoudre_une_instance("plat/exemple.dat")
# resoudre_une_instance("plat/plat10.dat")
# resoudre_une_instance("relief/relief10.dat")