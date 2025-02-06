import numpy as np
import pandas as pd

import plotly.graph_objects as go
import plotly.subplots
import plotly.io as pio
pio.renderers.default = 'png'

from typing import List, Union, Optional
from IPython.display import display, Math







# ====================================================================================================
# PHARMACOMETRICS ====================================================================================
# ====================================================================================================

# TRAPEIZOIDAL RULE ==================================================================================
def trapezoidal_rule(x, y):
    n = len(x)
    return sum((y[i] + y[i+1]) * (x[i+1] - x[i]) / 2 for i in range(n - 1))


# HENDERSON HASSELBALCH ==============================================================================
def henderson_hasselbalch(pKa=None, pH=None, fraction=None, type:str='ph', acid=True, printing=True):
    if (acid == True) and (printing == True):
        display(Math(r'\text{pH} = \text{pKa} + \log_{10} \left( \frac{\text{Ionized concentration}}{\text{Un-ionized concentration}} \right)'))
    elif (acid == False) and (printing == True):
        display(Math(r'\text{pH} = \text{pKa} + \log_{10} \left( \frac{\text{Un-ionized concentration}}{\text{Ionized concentration}} \right)'))
    if (type == 'ph'):
        return pKa+np.log10(fraction)
    if (type == 'pka'):
        return pH-np.log10(fraction)
    if (type == 'fraction'):
        return 10**(pH-pKa)
    


# RESIDUAL METHOD ====================================================================================
def residual_method(df:pd.DataFrame, elimination_end:int=5, absorption_start:int=5, show_fig:bool=True, parameters:bool=True) -> pd.DataFrame:

    fig = plotly.subplots.make_subplots(cols=2, rows=1,
        subplot_titles=['C x t','Residual Method']
    )

    # Conventional Plot
    fig.add_trace(trace=go.Scatter(
        x=df['t'], y=df['C'], mode='lines+markers', name='C x t',
        marker=dict(color='rgba(164,31,54,1)', size=8), line=dict(color='rgba(0,0,0,1)', shape='spline', smoothing=0.7)
    ), col=1, row=1)

    # Semi Log Plot
    fig.add_trace(trace=go.Scatter(
        x=df['t'], y=df['C'], mode='markers', name='semilog',
        marker=dict(color='rgba(0,0,0,1)', size=8),
    ), col=2, row=1)

    # Regression for Elimiation
    df['log10(C)'] = np.log10(df['C'])
    slope, intercept = np.polyfit(df['t'][-elimination_end:],df['log10(C)'][-elimination_end:],1)
    df['elimination line'] = 10**(df['t']*slope+intercept)
    fig.add_trace(trace=go.Scatter(
        x=df['t'], y=df['elimination line'], mode='lines', name='elimination',
        line=dict(color='rgba(164,31,54,1)'),
    ), col=2, row=1)

    # Residual Method
    df['residual'] = df['elimination line'] - df['C']
    absorption_slope, absorption_intercept = np.polyfit(df['t'][:absorption_start], np.log10(np.abs(df['residual']))[:absorption_start], 1)
    fig.add_trace(trace=go.Scatter(
        x=df['t'][:absorption_start], y=df['residual'][:absorption_start], mode='lines+markers', name='residual',
        line=dict(color='rgba(164,31,54,1)', dash='dash'), marker=dict(color='rgba(164,31,54,1)', size=8)
    ), col=2, row=1)

    if show_fig:
        main_subplot_layout(fig, title='Absorption Kinetics', x='t [h]', y='C [mg/L]', yaxis2_type="log").show()

    if parameters:
        k = round(-slope * np.log(10),2)
        half_life = round(np.log(2)/k,2)
        ka = round(-absorption_slope*np.log(10),2)
        absorption_half_life = round(np.log(2)/ka,2)
        AUC = round(trapezoidal_rule(x=df['t'], y=df['C']),2)
        return {'k':k, 'ka':ka, 't(1/2)':half_life, 't(1/2)(absorption)':absorption_half_life, 'AUC':AUC}

    return df



# WAGNER-NELSON METHOD ===============================================================================
def wagner_nelson(df:pd.DataFrame, half_life:float=None, k:float=None) -> pd.DataFrame:

    if (k is None):
        k = np.log(2)/half_life
    if (half_life is None):
        half_life = np.log(2)/half_life

    df['AUC'] = 0.
    for i in range(1, len(df)):
        df.loc[i, 'AUC'] = round(df.loc[i - 1, 'AUC'] + 0.5 * (df.loc[i, 'C'] + df.loc[i - 1, 'C']) * (df.loc[i, 't'] - df.loc[i - 1, 't']),2)
    df['A(ab)/V'] = round(df['C'] + k*df['AUC'],2)
    df['Bioavailable Drug Absorbed'] = round(df['A(ab)/V']/(k*df['AUC'][-1:].values),3)
    df['Bioavailable Drug Remained to be Absorbed'] = round(1-df['Bioavailable Drug Absorbed'],3)

    return df



# EQUATIONS ==========================================================================================
class Equations:
    def __init__(self):
        ...

    def chap_3(self, print_all=False, show_attr=False):
        self.extraction = Math(r'\text{E} = \frac{\left( \text{C} - \text{C}_{\text{out}} \right)}{\text{C}}')
        self.rate_of_elimination = Math(r'\text{Rate of elimination} = \text{CL} \cdot \text{C}')
        self.clearance_q_e = Math(r'\text{CL} = \text{Q} \cdot \text{E}')
        self.k_elimination_rate_constant = Math(r'\text{k} = \frac{\text{CL}}{\text{V}}')
        self.rate_of_elimination_differential_iv_bolus = Math(r'\text{Rate of elimination} = -\frac{d\text{A}}{dt} = \text{k} \cdot \text{A}')
        self.concentration_a_v = Math(r'\text{C} = \frac{\text{A}}{\text{V}}')
        self.concentration_time_elimination_iv_bolus = Math(r'\text{C} = \text{C}(0) \cdot e^{-\text{k} \cdot \text{t}}')
        self.log_concentration_time_elimination_iv_bolus = Math(r'\ln \text{C} = \ln \text{C}(0) - \text{k} \cdot \text{t}')
        self.half_life_k = Math(r'\text{t}_{1/2} = \frac{\ln 2}{\text{k}}')
        self.half_life_v_cl = Math(r'\text{t}_{1/2} = \frac{\ln 2 \cdot \text{V}}{\text{CL}}')
        self.fraction_dose_remaining = Math(r'\text{Fraction of dose remaining} = e^{-\text{k} \cdot \text{t}} = \left( \frac{1}{2} \right)^{\text{n}}')
        self.dose_cl_auc = Math(r'\text{Dose} = \text{CL} \cdot \text{AUC}')
        self.volume_distribution_a_c = Math(r'\text{V} = \frac{\text{A}}{\text{C}}')
        self.volume_distribution_half_life = Math(r'\text{V} = 1.44 \cdot \text{CL} \cdot \text{t}_{1/2}')
        self.mean_residence_time = Math(r'\text{MRT} = \frac{1}{\text{k}}')
        self.fe = Math(r'f_\text{e} = \frac{\text{A}_e\infty}{\text{Dose}}')
        self.renal_clearance_fe_cl = Math(r'\text{CL}_\text{R} = f\text{e} \cdot \text{CL}')
        self.renal_clearance_rate_excretion = Math(r'\text{CL}_\text{R} = \frac{\text{Rate of excretion}}{\text{Plasma concentration}}')
        self.renal_clearance_auc = Math(r'\text{CL}_\text{R} = \frac{\text{A}_e\infty}{\text{AUC}}')
        self.amount_excreted_iv_bolus = Math(r'\text{A}_e = \text{A}_e\infty \cdot e^{-\text{k} \cdot \text{t}}')
        self.fraction_eliminated = Math(r'\text{Fraction Eliminated} = \frac{\text{AUC}_{(0-x)}}{\text{AUC}_{(0-\infty)}}')
        if print_all == True:
            for attr, value in vars(self).items():
                if isinstance(value, Math):
                    if show_attr == True: display(attr, value); print()
                    else: display(value)

    def chap_4(self, print_all=False, show_attr=False):
        self.net_rate_transport = Math(r'\text{Net Rate of Transport} = \text{P} \cdot \text{SA} \cdot \left( \text{Cu}_{1} - \text{Cu}_{2} \right)')
        self.distribution_half_life = Math(r'\text{Distribution Half-Life} = \frac{\ln 2 \cdot \text{K}_{\text{p}_\text{b}}}{\text{Q}/\text{V}_\text{T}}')
        self.constant_equilibrium_distribution_rate = Math(r'\text{K}_{\text{p}_\text{b}} = \dfrac{\text{C}_\text{Tissue}}{\text{C}_\text{Blood}}')
        self.volume_c_a = Math(r'\text{V} = \frac{\text{A}}{\text{C}}')
        self.volume_compartiments = Math(r'\text{V} = \text{V}_\text{p} + \text{V}_{\text{T}_1} \cdot \text{K}_{\text{P}_1} + \text{V}_{\text{T}_2} \cdot \text{K}_{\text{P}_2}\> , \dots')
        self.fraction_inside_plasma = Math(r'\text{Fraction Inside Plasma} = \frac{\text{V}_\text{p}}{\text{V}}')
        self.fraction_outside_plasma = Math(r'\text{Fraction Outside Plasma} = \frac{\text{V} - \text{V}_\text{p}}{\text{V}}')
        self.amount_mass_balance = Math(r'\text{V} \cdot \text{C} = \text{V}\text{u} \cdot \text{Cu} = \text{V}_\text{b} \cdot \text{C}_\text{b}')
        self.fu = Math(r'f\text{u} = \frac{\text{Cu}}{\text{C}}')
        self.fu_protein = Math(r'f\text{u} = \frac{1}{1 + \text{Ka} \cdot f\text{u}_\text{p} \cdot \text{P}_\text{t}}')
        self.volume_tecidual_distribution = Math(r'\text{V} = \text{V}_\text{p} + \text{V}_\text{TW} \cdot \frac{f\text{u}}{f\text{u}_\text{T}}')
        self.partition_coefficient = Math(r'\text{K}_\text{p} = \frac{\text{C}_\text{T}}{\text{C}} = \frac{f\text{u}}{f\text{u}_\text{T}} \cdot \frac{\text{Cu}_\text{T}}{\text{Cu}}')
        self.volume_tecidual_distribution_permeability = Math(r'\text{V} = \text{V}_\text{p} + \text{V}_\text{TW} \cdot \frac{f\text{u}}{f\text{u}_\text{T}} \cdot \frac{\text{P}_\text{uptake}}{\text{P}_{\text{efflux}}}')
        self.volume_distribution_albumin = Math(r'\text{V} = 7.5 + \left( 7.5 + \frac{\text{V}_\text{R}}{f\text{u}_\text{R}} \right) \cdot f\text{u}')
        self.new_fu_estimation = Math(r"\dfrac{f\text{u}'}{f\text{u}} \approx \dfrac{\text{P}_\text{t}}{\text{P}_\text{t}'}")
        self.fu_r = Math(r"f\text{u}_\text{R} = \frac{\text{V}_\text{R}}{(\text{V} - 7.5) \cdot \frac{1}{f\text{u}}}")
        if print_all == True:
            for attr, value in vars(self).items():
                if isinstance(value, Math):
                    if show_attr == True: display(attr, value); print()
                    else: display(value)

    def chap_5(self, print_all=False, show_attr=False):
        self.renal_clearance_blood = Math(r'\text{CL}_{\text{R}} = \text{Q}_{\text{R}} \cdot \text{E}_{\text{R}}')
        self.renal_clearance_hepatic = Math(r'\text{CL}_{\text{H}} = \text{Q}_{\text{H}} \cdot \text{E}_{\text{H}}')
        self.rate_of_elimination = Math(r'\text{Rate of Elimination} = \text{CL} \cdot \text{C} = \text{CL}_\text{b} \cdot \text{C}_\text{b} = \text{CLu} \cdot \text{Cu}')
        self.clearance_b = Math(r'\text{CL}_{\text{b}} = \text{CL} \cdot \left(\frac{\text{C}}{\text{C}_{\text{b}}}\right)')
        self.clearance_u = Math(r'\text{CLu} = \frac{\text{CL}}{f\text{u}}')
        self.clearance = Math(r'\text{CL} = \text{CL}_{\text{R}} + \text{CL}_{\text{H}}')
        self.clearance_hepatic_blood = Math(r'\text{CL}_{\text{b,H}} = \text{Q}_{\text{H}} \cdot \left[\frac{f\text{u}_\text{b} \cdot \text{CL}_\text{int}}{\text{Q}_\text{H} + f\text{u}_\text{b} \cdot \text{CL}_\text{int}}\right]')
        self.extraction_ratio_hepatic = Math(r'\text{E}_{\text{H}} = \frac{f\text{u}_{\text{b}} \cdot \text{CL}_{\text{int}}}{\text{Q}_{\text{H}} + f\text{u}_{\text{b}} \cdot \text{CL}_{\text{int}}}')
        self.biliary_clearance = Math(r'\text{Biliary Clearance} = \frac{(\text{Bile Flow}) \cdot (\text{Concentration in Bile})}{\text{Concentration in Plasma}}')
        self.renal_clearance_rate = Math(r'\text{CL}_{\text{R}} = (1 - \text{F}_{\text{R}}) \cdot \left[\frac{\text{Rate of Filtration}}{\text{C}} + \frac{\text{Rate of Secretion}}{\text{C}}\right] = (1 - \text{F}_{\text{R}}) \cdot \left[\text{CL}_{\text{f}} + \text{CL}_{\text{s}}\right]')
        self.renal_clearance_excretion = Math(r'\text{CL}_{\text{R}} = \frac{\text{Rate of Urinary Excretion}}{\text{Plasma Concentration}}')
        self.renal_clearance_fu = Math(r'\text{Renal Clearance} = f\text{u} \cdot \text{Urine Flow}')
        self.k_cl_v = Math(r'\text{k} = \frac{\text{Rate of Elimination}}{\text{Amount in Body}} = \frac{\text{CL}}{\text{V}}')
        self.half_life = Math(r'\text{t}_{1/2} = \frac{\ln 2}{\text{k}} = \frac{\ln 2 \cdot \text{V}}{\text{CL}}')
        self.clearance_volume_relation = Math(r'\text{k} = \frac{\text{CL}}{\text{V}} = \frac{\text{CL}_{\text{b}}}{\text{V}_{\text{b}}} = \frac{\text{CLu}}{\text{Vu}}')
        self.half_life_clearance_relation = Math(r'\text{t}_{1/2} = \frac{\ln 2 \cdot \text{V}}{\text{CL}} = \frac{\ln 2 \cdot \text{V}_{\text{b}}}{\text{CL}_{\text{b}}} = \frac{\ln 2 \cdot \text{Vu}}{\text{CLu}}')
        if print_all == True:
            for attr, value in vars(self).items():
                if isinstance(value, Math):
                    if show_attr == True: display(attr, value); print()
                    else: display(value)

    def appendix_B(self, print_all=False, show_attr=False):
        self.henderson_hasselbalch_acid = Math(r'\text{pH} = \text{pKa} + \log_{10} \left( \frac{\text{Ionized concentration}}{\text{Un-ionized concentration}} \right)')
        self.henderson_hasselbalch_basic = Math(r'\text{pH} = \text{pKa} + \log_{10} \left( \frac{\text{Un-ionized concentration}}{\text{Ionized concentration}} \right)')
        if print_all == True:
            for attr, value in vars(self).items():
                if isinstance(value, Math):
                    if show_attr == True: display(attr, value); print()
                    else: display(value)
# ====================================================================================================





# ====================================================================================================
# PLOTINGS ===========================================================================================
# ====================================================================================================

# USED WITH PLOTLY TO ELABORATE THE LAYOUT ===========================================================
def main_layout(fig:go.Figure, width=700, height=600, x=None, y=None, title=None,
               x_range=None, y_range=None, paper_color='white', 
               customdata:Union[str, None]=None, hover_customdata='Info', 
               hover_x='x',hover_y='y', **kwargs) -> go.Figure:
    fig.layout = go.Layout(
        width=width,
        height=height,
        plot_bgcolor=paper_color,
        paper_bgcolor=paper_color,
        xaxis={'gridcolor':'#cccccc', 'linecolor':'black','title':x, 'range':x_range},
        yaxis={'gridcolor':'#cccccc', 'linecolor':'black','title':y, 'range':y_range},
        title={'text':title},
        **kwargs
    )
    if customdata == 'no':
        ...
    elif customdata is None:
        fig.update_traces(patch={
            'customdata':customdata, 'hovertemplate': hover_x + ': %{x}<br>' + hover_y + ': %{y}'
        })
    else:
        fig.update_traces(patch={
            'customdata':customdata,
            'hovertemplate': hover_x + ': %{x}<br>' + hover_y + ': %{y}<br>' + hover_customdata + ': %{customdata}<br>'
        })
    return fig
# ====================================================================================================
def main_subplot_layout(fig:go.Figure, width=1400, height=500, title=None, paper_color='white',
                        x=None, y=None, rows=1, cols=2, x_range=None, y_range=None,
                        customdata:Union[str, None]=None, hover_customdata='Info', 
                        hover_x='x',hover_y='y', **kwargs) -> go.Figure:
    fig.update_layout({
        'width':width,
        'height':height,
        'plot_bgcolor':paper_color,
        'paper_bgcolor':paper_color,
        'title':title,
        **kwargs
    })
    for xaxis in fig.select_xaxes():
        xaxis.update(
            showgrid=True,
            gridcolor='#CCCCCC',
            linecolor='black',
            title=x,
            range=x_range
        )
    for yaxis in fig.select_yaxes():
        yaxis.update(
            showgrid=True,
            gridcolor='#CCCCCC',
            linecolor='black',
            title=y,
            range=y_range
        )
    if customdata == 'no':
        ...
    elif customdata is None:
        fig.update_traces(patch={
            'customdata':customdata, 'hovertemplate': hover_x + ': %{x}<br>' + hover_y + ': %{y}'
        })
    else:
        fig.update_traces(patch={
            'customdata':customdata,
            'hovertemplate': hover_x + ': %{x}<br>' + hover_y + ': %{y}<br>' + hover_customdata + ': %{customdata}<br>'
        })
    return fig
# ====================================================================================================