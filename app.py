import numpy as np
import pandas as pd
import streamlit as st
from scipy.stats import linregress
import matplotlib.pyplot as plt

STP_ML_PER_MOL = 22414.0

st.set_page_config(page_title="NH3 GC Reaction Order Tool", layout="wide")
st.title("NH3 GC Reaction Order Calculator")

st.markdown("""
Reaction: **NH3 → 0.5 N2 + 1.5 H2**

Calibration equation: **Area = slope × x + intercept**  
where x is mole fraction or partial pressure. If x is partial pressure, set total pressure correctly.
""")

total_pressure = st.number_input("Total pressure for converting partial pressure to mole fraction (atm)", value=1.0, min_value=0.0001)
basis = st.selectbox("Calibration x-axis basis", ["mole fraction", "partial pressure"])

def default_cal(species):
    x = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
    factor = {"H2": 1000, "N2": 800, "NH3": 1200}[species]
    intercept = {"H2": 5, "N2": 3, "NH3": 8}[species]
    return pd.DataFrame({"x": x, "area": factor*x + intercept})

st.header("1. Calibration")
cols = st.columns(3)
cal_params = {}
for col, sp in zip(cols, ["H2", "N2", "NH3"]):
    with col:
        st.subheader(sp)
        df = st.data_editor(default_cal(sp), num_rows="dynamic", key=f"cal_{sp}")
        clean = df.dropna()
        if len(clean) >= 2:
            slope, intercept, r_value, _, _ = linregress(clean["x"], clean["area"])
            cal_params[sp] = {"slope": slope, "intercept": intercept, "r2": r_value**2}
            st.write(f"Area = {slope:.6g} × x + {intercept:.6g}; R²={r_value**2:.5f}")
        else:
            cal_params[sp] = None
            st.warning("Need at least two calibration points.")

st.header("2. Experiments")
default_exp = pd.DataFrame({
    "Run": [1,2,3,4,5],
    "Order target": ["H2"]*5,
    "Wcat_g": [0.1]*5,
    "Ftot_in_sccm": [100]*5,
    "F_H2_in": [5,10,15,20,25],
    "F_N2_in": [20]*5,
    "F_NH3_in": [25]*5,
    "F_Ar_in": [50,45,40,35,30],
    "Area_H2": [np.nan]*5,
    "Area_N2": [np.nan]*5,
    "Area_NH3": [np.nan]*5,
})
exp = st.data_editor(default_exp, num_rows="dynamic", key="exp")

def area_to_y(area, sp):
    p = cal_params.get(sp)
    if p is None or pd.isna(area):
        return np.nan
    x_axis = (area - p["intercept"]) / p["slope"]
    return x_axis if basis == "mole fraction" else x_axis / total_pressure

def calculate(row):
    y_h2 = area_to_y(row["Area_H2"], "H2")
    y_n2 = area_to_y(row["Area_N2"], "N2")
    y_nh3 = area_to_y(row["Area_NH3"], "NH3")
    ft = row["Ftot_in_sccm"]
    fh2, fn2, fnh3 = row["F_H2_in"], row["F_N2_in"], row["F_NH3_in"]

    x_h2 = (fh2 - y_h2*ft)/(y_h2 - 1.5) if pd.notna(y_h2) else np.nan
    x_n2 = (fn2 - y_n2*ft)/(y_n2 - 0.5) if pd.notna(y_n2) else np.nan
    x_nh3 = (fnh3 - y_nh3*ft)/(1 + y_nh3) if pd.notna(y_nh3) else np.nan
    x_avg = np.nanmean([x_h2, x_n2, x_nh3]) if not np.all(pd.isna([x_h2, x_n2, x_nh3])) else np.nan

    rate = x_avg / (STP_ML_PER_MOL * 60 * row["Wcat_g"]) if pd.notna(x_avg) and row["Wcat_g"] > 0 else np.nan
    p_h2 = fh2/ft*total_pressure
    p_n2 = fn2/ft*total_pressure
    p_nh3 = fnh3/ft*total_pressure
    target = row["Order target"]
    p_target = {"H2": p_h2, "N2": p_n2, "NH3": p_nh3}.get(target, np.nan)

    return pd.Series({
        "y_H2": y_h2, "y_N2": y_n2, "y_NH3": y_nh3,
        "x_from_H2_sccm": x_h2, "X_from_H2": x_h2/fnh3 if pd.notna(x_h2) else np.nan,
        "x_from_N2_sccm": x_n2, "X_from_N2": x_n2/fnh3 if pd.notna(x_n2) else np.nan,
        "x_from_NH3_sccm": x_nh3, "X_from_NH3": x_nh3/fnh3 if pd.notna(x_nh3) else np.nan,
        "x_avg_sccm": x_avg, "X_avg": x_avg/fnh3 if pd.notna(x_avg) else np.nan,
        "rate_mol_g_s": rate,
        "P_H2_in_atm": p_h2, "P_N2_in_atm": p_n2, "P_NH3_in_atm": p_nh3,
        "ln_rate": np.log(rate) if pd.notna(rate) and rate > 0 else np.nan,
        "ln_P_target": np.log(p_target) if pd.notna(p_target) and p_target > 0 else np.nan,
    })

if all(cal_params.values()):
    result = pd.concat([exp, exp.apply(calculate, axis=1)], axis=1)
    st.header("3. Results")
    st.dataframe(result, use_container_width=True)

    st.header("4. Reaction order")
    summary = []
    for target in ["H2", "N2", "NH3"]:
        sub = result[(result["Order target"] == target) & result["ln_rate"].notna() & result["ln_P_target"].notna()]
        if len(sub) >= 2:
            slope, intercept, r_value, _, _ = linregress(sub["ln_P_target"], sub["ln_rate"])
            summary.append({"Target": target, "Reaction order": slope, "ln(k_app)": intercept, "R2": r_value**2, "n": len(sub)})
            fig, ax = plt.subplots()
            ax.scatter(sub["ln_P_target"], sub["ln_rate"])
            xline = np.linspace(sub["ln_P_target"].min(), sub["ln_P_target"].max(), 50)
            ax.plot(xline, slope*xline + intercept)
            ax.set_xlabel(f"ln(P_{target})")
            ax.set_ylabel("ln(rate)")
            ax.set_title(f"{target} order = {slope:.3g}")
            st.pyplot(fig)
        else:
            summary.append({"Target": target, "Reaction order": np.nan, "ln(k_app)": np.nan, "R2": np.nan, "n": len(sub)})
    st.dataframe(pd.DataFrame(summary), use_container_width=True)

    csv = result.to_csv(index=False).encode("utf-8-sig")
    st.download_button("Download result CSV", csv, "nh3_gc_results.csv", "text/csv")
else:
    st.info("Enter valid calibration data for H2, N2, and NH3.")
