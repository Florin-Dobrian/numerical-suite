import numpy as np
import pandas as pd
import doubleml as dml
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor

# --- Nonlinear Confounding Experiment ---
# TRUE TREATMENT EFFECT = 3.0

np.random.seed(42)
n = 10000

# Generate 5 confounders
X1 = np.random.normal(0, 1, n)
X2 = np.random.normal(0, 1, n)
X3 = np.random.normal(0, 1, n)
X4 = np.random.normal(0, 1, n)
X5 = np.random.normal(0, 1, n)

# Treatment assignment: nonlinear function of confounders
treatment_propensity = (
    0.5 * X1**2 
    + 0.8 * np.sin(X2 * np.pi) 
    + 0.6 * X3 * X4
    + 0.3 * np.abs(X5)
    + np.random.normal(0, 0.5, n)
)
prob_treat = 1 / (1 + np.exp(-treatment_propensity))
D = np.random.binomial(1, prob_treat)

# Outcome: TRUE EFFECT of D is 3.0
Y = (
    3.0 * D
    + 2.0 * X1**2
    + 3.0 * np.sin(X2 * np.pi)
    + 2.5 * X3 * X4
    + 1.5 * np.abs(X5)
    + 1.0 * X1 * X2
    + np.random.normal(0, 1, n)
)

df = pd.DataFrame({
    'Y': Y, 'D': D,
    'X1': X1, 'X2': X2, 'X3': X3, 'X4': X4, 'X5': X5
})

x_cols = ['X1', 'X2', 'X3', 'X4', 'X5']
dml_data = dml.DoubleMLData(df, 'Y', 'D', x_cols=x_cols)

print("=" * 70)
print("Nonlinear Confounding Experiment (n=10,000)")
print("TRUE TREATMENT EFFECT = 3.0")
print("=" * 70)

# --- Elastic Net ---
X_en = df[['D', 'X1', 'X2', 'X3', 'X4', 'X5']]
en_model = ElasticNetCV(cv=5).fit(X_en, df['Y'])
en_effect = en_model.coef_[0]

# --- DoubleML with Random Forest ---
ml_l_rf = RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=5)
ml_m_rf = RandomForestRegressor(n_estimators=500, max_depth=None, min_samples_leaf=5)
dml_rf = dml.DoubleMLPLR(dml_data, ml_l_rf, ml_m_rf, n_folds=5)
dml_rf.fit()

# --- DoubleML with XGBoost ---
ml_l_xgb = XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.1, verbosity=0)
ml_m_xgb = XGBRegressor(n_estimators=500, max_depth=6, learning_rate=0.1, verbosity=0)
dml_xgb = dml.DoubleMLPLR(dml_data, ml_l_xgb, ml_m_xgb, n_folds=5)
dml_xgb.fit()

# --- DoubleML with LightGBM ---
ml_l_lgb = LGBMRegressor(n_estimators=500, max_depth=-1, learning_rate=0.1, verbose=-1)
ml_m_lgb = LGBMRegressor(n_estimators=500, max_depth=-1, learning_rate=0.1, verbose=-1)
dml_lgb = dml.DoubleMLPLR(dml_data, ml_l_lgb, ml_m_lgb, n_folds=5)
dml_lgb.fit()

# --- Summary ---
results = [
    ("Ground Truth",            3.000, 0.000, "-"),
    ("Elastic Net",             en_effect, en_effect - 3.0, "-"),
    ("DoubleML (Random Forest)", dml_rf.coef[0], dml_rf.coef[0] - 3.0,
     f"[{dml_rf.confint().iloc[0,0]:.3f}, {dml_rf.confint().iloc[0,1]:.3f}]"),
    ("DoubleML (XGBoost)",      dml_xgb.coef[0], dml_xgb.coef[0] - 3.0,
     f"[{dml_xgb.confint().iloc[0,0]:.3f}, {dml_xgb.confint().iloc[0,1]:.3f}]"),
    ("DoubleML (LightGBM)",     dml_lgb.coef[0], dml_lgb.coef[0] - 3.0,
     f"[{dml_lgb.confint().iloc[0,0]:.3f}, {dml_lgb.confint().iloc[0,1]:.3f}]"),
]

print(f"\n{'Method':<30} {'Estimate':>10} {'Bias':>10}   {'95% CI'}")
print("-" * 70)
for name, est, bias, ci in results:
    print(f"{name:<30} {est:>10.3f} {bias:>+10.3f}   {ci}")
print("=" * 70)
