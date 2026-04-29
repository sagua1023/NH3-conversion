# NH3 GC Reaction Order Tool

GC area data를 calibration curve로 mole fraction/partial pressure로 변환하고,
NH3 분해 반응의 전환율과 reaction order를 계산하는 간단한 Streamlit 앱입니다.

## Reaction
NH3 -> 0.5 N2 + 1.5 H2

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Inputs
1. Calibration data: H2, N2, NH3의 calibration x 값과 GC area
2. Experiment data: inlet flow, measured GC area, catalyst mass
3. Order target: H2/N2/NH3 중 어떤 분압을 변화시킨 실험인지 선택

## Core equations
- H2 기준: x = (F_H2,in - y_H2 F_tot,in)/(y_H2 - 1.5)
- N2 기준: x = (F_N2,in - y_N2 F_tot,in)/(y_N2 - 0.5)
- NH3 기준: x = (F_NH3,in - y_NH3 F_tot,in)/(1 + y_NH3)
- Conversion: X_NH3 = x/F_NH3,in
- Rate: r = x/(22414*60*Wcat)
- Reaction order: slope of ln(rate) vs ln(P_target)
