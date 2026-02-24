import numpy as np
import pandas as pd
import doubleml as dml
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor

# --- Nonlinear Confounding Experiment ---
# TRUE TREATMENT EFFECT = 3.0
# 
# The twist: confounders affect outcome and treatment through
# nonlinear functions (interactions, sine, squared terms).
# Elastic Net (linear) can't capture these relationships,
# so it produces a biased estimate. DoubleML with Random Forests can.

np.random.seed(42)
n = 5000

# Generate 5 confounders
X1 = np.random.normal(0, 1, n)
X2 = np.random.normal(0, 1, n)
X3 = np.random.normal(0, 1, n)
X4 = np.random.normal(0, 1, n)
X5 = np.random.normal(0, 1, n)

# Treatment assignment: nonlinear function of confounders
# Includes interactions and nonlinear transformations
treatment_propensity = (
    0.5 * X1**2 
    + 0.8 * np.sin(X2 * np.pi) 
    + 0.6 * X3 * X4  # interaction
    + 0.3 * np.abs(X5)
    + np.random.normal(0, 0.5, n)
)
prob_treat = 1 / (1 + np.exp(-treatment_propensity))
D = np.random.binomial(1, prob_treat)

# Outcome: TRUE EFFECT of D is 3.0
# Confounders affect outcome nonlinearly too
Y = (
    3.0 * D                          # true causal effect
    + 2.0 * X1**2                    # quadratic
    + 3.0 * np.sin(X2 * np.pi)      # sinusoidal
    + 2.5 * X3 * X4                  # interaction
    + 1.5 * np.abs(X5)              # absolute value
    + 1.0 * X1 * X2                  # another interaction
    + np.random.normal(0, 1, n)      # noise
)

df = pd.DataFrame({
    'Y': Y, 'D': D,
    'X1': X1, 'X2': X2, 'X3': X3, 'X4': X4, 'X5': X5
})

print("=" * 60)
print("Nonlinear Confounding Experiment")
print("TRUE TREATMENT EFFECT = 3.0")
print("=" * 60)

# --- Elastic Net (linear model) ---
X_en = df[['D', 'X1', 'X2', 'X3', 'X4', 'X5']]
y_en = df['Y']
en_model = ElasticNetCV(cv=5).fit(X_en, y_en)
print(f"\nElastic Net Estimated Effect:  {en_model.coef_[0]:.3f}")
print(f"  (Bias: {en_model.coef_[0] - 3.0:+.3f})")

# --- DoubleML with Random Forests ---
x_cols = ['X1', 'X2', 'X3', 'X4', 'X5']
dml_data = dml.DoubleMLData(df, 'Y', 'D', x_cols=x_cols)

ml_l = RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=5)
ml_m = RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=5)

dml_plr = dml.DoubleMLPLR(dml_data, ml_l, ml_m, n_folds=5)
dml_plr.fit()
print(f"\nDoubleML Estimated Effect:     {dml_plr.coef[0]:.3f}")
print(f"  (Bias: {dml_plr.coef[0] - 3.0:+.3f})")
print(f"  95% CI: [{dml_plr.confint().iloc[0, 0]:.3f}, {dml_plr.confint().iloc[0, 1]:.3f}]")

# --- Summary ---
print(f"\n{'=' * 60}")
print(f"{'Method':<30} {'Estimate':>10} {'Bias':>10}")
print(f"{'-' * 60}")
print(f"{'Ground Truth':<30} {'3.000':>10} {'0.000':>10}")
print(f"{'Elastic Net':<30} {en_model.coef_[0]:>10.3f} {en_model.coef_[0] - 3.0:>+10.3f}")
print(f"{'DoubleML (Random Forest)':<30} {dml_plr.coef[0]:>10.3f} {dml_plr.coef[0] - 3.0:>+10.3f}")
print(f"{'=' * 60}")
