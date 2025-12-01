
function clean_names!(d)
    rename!(lowercase, d)
    rename!(d, names(d) .=> replace.(names(d), r"\s+" => "_"))
    rename!(d, names(d) .=> replace.(names(d), r"," => ""))
    rename!(d, names(d) .=> replace.(names(d), r"\+" => ""))
    rename!(d, names(d) .=> replace.(names(d), r"\." => "_"))
    rename!(d, names(d) .=> replace.(names(d), r"__" => "_"))
    rename!(d, names(d) .=> replace.(names(d), r"\?" => ""))
end


function select_samples!(d)
    filter!(row -> row.vendor != "ICR", d)
    #filter!(row -> row.collection_country != "Unknown", d)
    #filter!(row -> row.stage in ["1", "2", "N/A"], d)
    #filter!(row -> row.race in ["White", "Hispanic"], d)
    filter!(row -> row.sex == "M", d)
    filter!(row -> row.disease in ["Healthy", "Prostate Cancer"], d)
    filter!(row -> row.training_or_blinded != "Z_Do_Not_Analyze_Hemolyzed", d)
    filter!(row -> !ismissing(row.total_psa_access), d)
    filter!(:serum_specific_id => !in(["HMN686019S", "HMN456748"]), d)
end



function select_variables!(d)
    select!(d, :serum_specific_id, :vendor, :targeted_disease_for_collection,
            :collection_country, :age_at_collection, :sex, :race, :stage,
            :total_psa_access, :free_psa_access, :p2psa_access, :phi_score_access,
            :cancer, :disease)

    rename!(d, [:serum_specific_id => :id,
                :targeted_disease_for_collection => :targeted_disease,
                :collection_country => :country,
                :age_at_collection => :age]
            )
    
    rename!(d, names(d) .=> replace.(names(d), r"_access" => ""))

    
    d.bound_psa = d.total_psa .- d.free_psa
    d.bound_psa[d.bound_psa .< 0.0] .= 0.0
    d.phi_score = tofloat.(d.phi_score)
    d.age = tofloat.(d.age)

    filter!(row -> 52.0 < row.age < 83.0, d)

end

jitter(x, factor = 0.03) = x + factor * randn(size(x))

function invlogit(x)
  1 / (1 + exp(-x))
end


function tofloat(str)
    tryparse(Float64, str) !== nothing ? parse(Float64, str) : missing
end



function extract_returned(model, chains::Chains)

    ret = returned(model, chains)

    N_iter, N_chains = size(ret)
    N_samples = length(ret[1, 1])

    tmp = vec(ret)

    mat = zeros(N_iter * N_chains, N_samples)
    for i in 1:(N_iter * N_chains)
        mat[i, :] = tmp[i]
    end
    return mat
end



function plot_pred(x, y; title = "")
    scatter(jitter(x), y, label = :none, title = title, 
            xlabel = "\nObserved Cancer Status",
            ylabel = "Predicted Probability of Cancer\n", 
            ylim = (-0.03, 1.03), xlim = (-0.5, 1.5), 
            xticks = ([0, 1], ["Healthy", "Cancer"]),
            guidefont = font(12),  # axis labels
            tickfont  = font(12),  # tick labels
            tick_direction = :out,
            markercolour = :lightgrey,
            markerstrokecolour = :steelblue,
            size = (500, 500))
end 
