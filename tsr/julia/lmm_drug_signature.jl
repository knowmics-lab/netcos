using DataFrames, MixedModels, StatsModels, CSV, Dates, LinearAlgebra, JLD2

tot_threads = Threads.nthreads()
experiments_meta_data_dir = "output/drug_signature/LINCS/experiments_data/"
signature_dir = "output/drug_signature/LINCS/julia_genes/"
gene_list_filename = "data/LINCS-GSE92742/julia_genes_to_compute.csv"

function map_drug_names(fixed_effect_names)
    fixed_effect_names_size = length(fixed_effect_names)
    drugs = Vector{String}(undef, fixed_effect_names_size - 1)
    for i in 2:fixed_effect_names_size
        drugs[i - 1] = chop(fixed_effect_names[i], head = 12, tail = 0)
    end
    return drugs
end

function map(LMM_output::LinearMixedModel{Float64}, gene_symbol::String)
    coef_size = length(LMM_output.beta)
    drug_signature = DataFrame(
        drugs = map_drug_names(LMM_output.feterm.cnames),
        gene = fill(gene_symbol,coef_size - 1),
        DE_log2_FC = LMM_output.beta[2:coef_size],
        std_error = LMM_output.stderror[2:coef_size],
        t_value = Vector{Float64}(undef,coef_size - 1),
        p_value = LMM_output.pvalues[2:coef_size]
    )
    drug_signature.t_value = drug_signature.DE_log2_FC./drug_signature.std_error
    return drug_signature
end

function lmm_compute(experiments_data_filename::String, gene_symbol::String, computation_number::Integer)
    #=
    if computation_number <= tot_threads
        println(string("sleep computation number: ",computation_number))
        sleep(10 * (computation_number - 1))
    end
    =#
    println(string("start compuation: ",computation_number,", filename: ",experiments_data_filename,", time: ",now()))
    experiments_data = CSV.read(experiments_meta_data_dir*experiments_data_filename*".csv",DataFrame)
    println(string("start compuation after read csv: ",computation_number,", time: ",now()))
    start_c = now()
    LMM_output = fit!(LinearMixedModel(@formula(gene_expression ~ pert_iname + (1 | cell_id) + (1 | rna_plate)), experiments_data), REML=true);
    end_c = now()
    println(string("end computation lmm: ",computation_number,", time: ",now(),", computaion time: ",end_c-start_c))
    LMM_output=map(LMM_output,gene_symbol)
    println(string("end compuation after map: ",computation_number,", time: ",now()))
    filename=signature_dir*experiments_data_filename*".jld2"
    println(filename)
    save_object(filename, LMM_output)
    return LMM_output
end

#BLAS.set_num_threads(4)

genes = CSV.read(gene_list_filename, DataFrame)

start_c = now()
println(string("start overall compuation time: ", now()))

tot_genes = length(genes.gene)
Threads.@threads for i = 1:tot_genes
   lmm_compute(String(genes.filename[i]),String(genes.gene[i]),i)
end

end_c = now()
println(string("end overall computation time: ",now(), ", total time: ",end_c - start_c))
