using Pkg
Pkg.activate(".")

using Turing
using Distributions
using StatsPlots
using Plots: PlotMeasures.px
using Random
using LinearAlgebra
using DataFrames
using CSV
using Splines2
using RCall
using MLUtils: splitobs
using GLM
using StatsBase
#using AnovaGLM
#using StatsFuns: logistic
using StatisticalMeasures

include("functions.jl")


d = CSV.read("../data/Cancer Model Training Data 20250828.csv", DataFrame)

clean_names!(d)
select_samples!(d)
select_variables!(d)

select(d, r"psa")


# println.(names(d))
# describe(d)
# CSV.write("../data/cleaned_data.csv", d)

# scatter(log10.(d.total_psa .+ 1), log10.(d.free_psa .+ 1), group=d.cancer)
# scatter(log10.(d.bound_psa .+ 1), log10.(d.free_psa .+ 1), group=d.cancer)
# scatter(log10.(d.total_psa .+ 1), log10.(d.p2psa .+ 1), group=d.cancer)

# boxplot(d.vendor, log10.(d.total_psa .+ 1), group=d.cancer)
# boxplot(d.vendor, log10.(d.free_psa .+ 1), group=d.cancer)
# boxplot(d.vendor, log10.(d.p2psa .+ 1), group=d.cancer)

# boxplot(d.country, log10.(d.total_psa .+ 1), group=d.cancer)
# scatter(d.country, log10.(d.total_psa .+ 1), group=d.cancer)
# boxplot(d.country, log10.(d.free_psa .+ 1), group=d.cancer)
# boxplot(d.country, log10.(d.p2psa .+ 1), group=d.cancer)

# boxplot(d.race, log10.(d.total_psa .+ 1), group=d.cancer)
# boxplot(d.race, log10.(d.free_psa .+ 1), group=d.cancer)
# boxplot(d.race, log10.(d.p2psa .+ 1), group=d.cancer)

# scatter(d.stage, log10.(d.total_psa .+ 1), group=d.cancer)
# scatter(d.stage, log10.(d.free_psa .+ 1), group=d.cancer)
# scatter(d.stage, log10.(d.p2psa .+ 1), group=d.cancer)


histogram(log10.(d.total_psa .+ 1))
histogram(log10.(d.bound_psa .+ 1))
histogram(log10.(d.free_psa .+ 1))
histogram(log10.(d.p2psa .+ 1))

## do train/test split
Random.seed!(123)
train_df, test_df = splitobs(d, at = 0.8; stratified = d.cancer)
train_df = DataFrame(train_df)
test_df =  DataFrame(test_df)




bound_norm_obj = R"bestNormalize::orderNorm($(train_df.bound_psa[train_df.cancer.==0]))"
free_norm_obj = R"bestNormalize::orderNorm($(train_df.free_psa[train_df.cancer.==0]))"
p2psa_norm_obj = R"bestNormalize::orderNorm($(train_df.p2psa[train_df.cancer.==0]))"

train_norm = hcat(rcopy(R"predict($bound_norm_obj, $(train_df.bound_psa))"), 
                  rcopy(R"predict($free_norm_obj, $(train_df.free_psa))"), 
                  rcopy(R"predict($p2psa_norm_obj, $(train_df.p2psa))")
                  )


test_norm = hcat(rcopy(R"predict($bound_norm_obj, $(test_df.bound_psa))"), 
                 rcopy(R"predict($free_norm_obj, $(test_df.free_psa))"), 
                 rcopy(R"predict($p2psa_norm_obj, $(test_df.p2psa))")
                 )


train_norm = DataFrame(train_norm, [:bound_psa, :free_psa, :p2psa])
train_norm.cancer = train_df.cancer
test_norm = DataFrame(test_norm, [:bound_psa, :free_psa, :p2psa])
test_norm.cancer = test_df.cancer


m1 = glm(@formula(cancer ~ bound_psa + free_psa + p2psa),
         train_norm, Binomial())

m1b = glm(@formula(cancer ~ bound_psa + free_psa * p2psa),
         train_norm, Binomial())
aicc(m1)
aicc(m1b)

cornerplot(hcat(train_norm.bound_psa, train_norm.free_psa, train_norm.p2psa))

mx0 = glm(hcat(train_norm.bound_psa), train_norm.cancer, Binomial())
#mx1 = glm(B1, train_norm.cancer, Binomial())
mx1 = glm(is(train_norm.bound_psa, df = 3, order = 3), train_norm.cancer, Binomial())
mx2 = glm(is(train_norm.bound_psa, df = 4, order = 3), train_norm.cancer, Binomial())

scatter(train_norm.bound_psa, logit.(predict(mx0)), group = train_norm.cancer)
scatter(train_norm.bound_psa, logit.(predict(mx1)), group = train_norm.cancer)
scatter(train_norm.bound_psa, logit.(predict(mx2)), group = train_norm.cancer)

aicc(mx0)
aicc(mx1) # best for bound psa
aicc(mx2)

mx0 = glm(hcat(train_norm.free_psa), train_norm.cancer, Binomial())
mx1 = glm(bs(train_norm.free_psa, df = 3, order = 3), train_norm.cancer, Binomial())
mx2 = glm(bs(train_norm.free_psa, df = 4, order = 3), train_norm.cancer, Binomial())

aicc(mx0)
aicc(mx1)
aicc(mx2) # best for free psa

mx0 = glm(hcat(train_norm.p2psa), train_norm.cancer, Binomial())
mx1 = glm(bs(train_norm.p2psa, df = 3, order = 3), train_norm.cancer, Binomial())
mx2 = glm(bs(train_norm.p2psa, df = 4, order = 3), train_norm.cancer, Binomial())

aicc(mx0)
aicc(mx1)
aicc(mx2) # best for p2psa




scatter(jitter(train_df.cancer), predict(m1), label = :none,
        xlabel = "Observed Cancer Status",
        ylabel = "Predicted Probability of Cancer"
        )

m1_test = predict(m1, test_norm)

scatter(jitter(test_df.cancer), m1_test, label = :none,
        xlabel = "Observed Cancer Status",
        ylabel = "Predicted Probability of Cancer", 
        ylim = (0, 1))


sum(train_df.cancer)
sum(test_df.cancer)

B1_obj = R"splines2::iSpline($(train_norm.bound_psa), knots=c(-1.0, 1.0), degree=3)"
B2_obj = R"splines2::iSpline($(train_norm.free_psa), knots=c(-1.0, 1.0), degree=3)"
B3_obj = R"splines2::iSpline($(train_norm.p2psa), knots=c(-1.0, 1.0), degree=3)"

# B1 = is(train_norm.bound_psa, interior_knots = [-1.0, 1.0], order = 3)
# B2 = is(train_norm.free_psa, interior_knots = [-1.0, 1.0], order = 3)
# B3 = is(train_norm.p2psa, interior_knots = [-1.0, 1.0], order = 3)

B1 = rcopy(B1_obj)
B2 = rcopy(B2_obj)
B3 = rcopy(B3_obj)

B1_test = rcopy(R"predict($B1_obj, $(test_norm.bound_psa))")
B2_test = rcopy(R"predict($B2_obj, $(test_norm.free_psa))")
B3_test = rcopy(R"predict($B3_obj, $(test_norm.p2psa))")

d_new = DataFrame(hcat(B1, B2, B3), :auto)
d_new.cancer = train_df.cancer

m2 = glm(@formula(cancer ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9),
         d_new, Binomial())


scatter(jitter(train_df.cancer), predict(m2), label = :none,
        xlabel = "Observed Cancer Status",
        ylabel = "Predicted Probability of Cancer"
        )

@model function prostate_mod(X₁, X₂, X₃)

    N, P = size(X₁)
    
    # priors
    β₀ ~ Normal(0, 50);

    σ₁_penalty ~ Exponential(0.5)
    σ₂_penalty ~ Exponential(0.5)
    σ₃_penalty ~ Exponential(0.5)

    β₁ ~ filldist(truncated(Normal(0, σ₁_penalty), 0, 100), P)
    β₂ ~ filldist(truncated(Normal(0, σ₂_penalty), 0, 100), P)
    β₃ ~ filldist(truncated(Normal(0, σ₃_penalty), 0, 100), P)

    μ = β₀ .+ X₁ * β₁ .+ X₂ * β₂ .+ X₃ * β₃
    y ~ product_distribution(BernoulliLogit.(μ))

    return μ
end
    
# model = prostate_mod(train_norm.cancer, B1, B2, B3)
# chains = sample(model, NUTS(), MCMCThreads(), 1_000, 3)
# describe(chains)

model = prostate_mod(B1, B2, B3)
model_cond = model | (y = train_norm.cancer, )
chains = sample(model_cond, NUTS(), MCMCThreads(), 1_000, 3)

chains_df = DataFrame(chains)

boxplot(Matrix(select(chains_df, r"σ")))

describe(chains)

pred_train = predict(model, chains)

mu_pred_train = extract_returned(model, chains)

pred = mean(invlogit.(mu_pred_train), dims = 1)[:]


savefig(
    plot_pred(train_norm.cancer, pred, title = "GAM model"),
    "figs/testplot.pdf"
)


model_test = prostate_mod(B1_test, B2_test, B3_test)
pred_test = predict(model_test, chains)
pred_test = extract_returned(model_test, chains)
pred_test = mean(invlogit.(pred_test), dims = 1)[:]

plot_pred(test_norm.cancer, pred_test, title = "GAM model")

describe(train_norm)

xs_bound = -1.5:0.05:5
xs_free = -2.4:0.05:2.8
xs_p2 = -2.8:0.05:2.8


## ============================================================
## Bound PSA
B1_plot = rcopy(R"predict($B1_obj, $(xs_bound))")
B2_plot = rcopy(R"predict($B2_obj, $(zeros(length(xs_bound))))")
B3_plot = rcopy(R"predict($B3_obj, $(zeros(length(xs_bound))))")

model_plot = prostate_mod(B1_plot, B2_plot, B3_plot)
pred_plot = extract_returned(model_plot, chains)

means = mean(pred_plot, dims = 1)[:]
cis = mapslices(col -> quantile(col, [0.025, 0.975]), pred_plot, dims = 1)

p1 = plot(xs_bound, means, label = :none,
     #ribbon = (cis[1, :], cis[2, :]), 
     xlabel = "\nBound PSA", 
     ylabel = "Log-odds of Cancer\n",
     guidefont = font(12),  # axis labels
     tickfont  = font(12),  # tick labels
     tick_direction = :out,
     linecolour = :steelblue,
     linewidth = 2
     #size = (500, 500)
     )
#plot!(xs, cis[1, :])
#plot!(xs, cis[2, :])

scatter!(train_norm.bound_psa[train_norm.cancer .== 0], fill(minimum(means) - 0.1, 166),
         label = :none, markershape = :vline, color = :black)

scatter!(train_norm.bound_psa[train_norm.cancer .== 1], fill(maximum(means) + 0.1, 41),
         label = :none, markershape = :vline, color = :black)


## ============================================================
## Free PSA

B1_plot = rcopy(R"predict($B1_obj, $(zeros(length(xs_free))))")
B2_plot = rcopy(R"predict($B2_obj, $(xs_free))")
B3_plot = rcopy(R"predict($B3_obj, $(zeros(length(xs_free))))")


model_plot = prostate_mod(B1_plot, B2_plot, B3_plot)
pred_plot = extract_returned(model_plot, chains)

means = mean(pred_plot, dims = 1)[:]
cis = mapslices(col -> quantile(col, [0.025, 0.975]), pred_plot, dims = 1)

p2 = plot(xs_free, means, label = :none,
     #ribbon = (cis[1, :], cis[2, :]), 
     xlabel = "\nFree PSA", 
     ylabel = "Log-odds of Cancer\n",
     guidefont = font(12),  # axis labels
     tickfont  = font(12),  # tick labels
     tick_direction = :out,
     linecolour = :steelblue,
     linewidth = 2
     #size = (500, 500)
     )
#plot!(xs, cis[1, :])
#plot!(xs, cis[2, :])

scatter!(train_norm.free_psa[train_norm.cancer .== 0], fill(minimum(means) - 0.1, 166),
         label = :none, markershape = :vline, color = :black)

scatter!(train_norm.free_psa[train_norm.cancer .== 1], fill(maximum(means) + 0.1, 41),
         label = :none, markershape = :vline, color = :black)


## ============================================================
## P2PSA

B1_plot = rcopy(R"predict($B1_obj, $(zeros(length(xs_p2))))")
B2_plot = rcopy(R"predict($B2_obj, $(zeros(length(xs_p2))))")
B3_plot = rcopy(R"predict($B3_obj, $(xs_p2))")


model_plot = prostate_mod(B1_plot, B2_plot, B3_plot)
pred_plot = extract_returned(model_plot, chains)

means = mean(pred_plot, dims = 1)[:]
cis = mapslices(col -> quantile(col, [0.025, 0.975]), pred_plot, dims = 1)

p3 = plot(xs_p2, means, label = :none,
     #ribbon = (cis[1, :], cis[2, :]), 
     xlabel = "\nP2PSA", 
     ylabel = "Log-odds of Cancer\n",
     guidefont = font(12),  # axis labels
     tickfont  = font(12),  # tick labels
     tick_direction = :out,
     linecolour = :steelblue,
          linewidth = 2
     #size = (500, 500)
     )
#plot!(xs, cis[1, :])
#plot!(xs, cis[2, :])

scatter!(train_norm.p2psa[train_norm.cancer .== 0], fill(minimum(means) - 0.1, 166),
         label = :none, markershape = :vline, color = :black)

scatter!(train_norm.p2psa[train_norm.cancer .== 1], fill(maximum(means) + 0.1, 41),
         label = :none, markershape = :vline, color = :black)


savefig(
    plot(p1, p2, p3, layout = (1, 3),
     left_margin = 40px, 
     bottom_margin = 40px, 
     windowsize = (400 * 3, 400) # width, height
     ), 
    "figs/prostate_preds.pdf"
)


function bs(ŷ, y)
    sum((ŷ .- y).^2)
end

bs(pred, train_norm.cancer)

StatisticalMeasures.Functions.auc(pred, train_norm.cancer, 1)
BalancedAccuracy()(pred .> 0.5, train_norm.cancer)
Accuracy()(pred .> 0.5, train_norm.cancer)
cm = confmat(pred .> 0.5, train_norm.cancer)

sensitivity(cm)
specificity(cm)
ppv(cm)
npv(cm)

## HERE, need to add calibration metrics
