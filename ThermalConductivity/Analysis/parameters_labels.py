# Creation of a dictionnary to sort data
def get_parameters_for_measures():
    dict_measures = dict()
    dict_measures["T0"] = ["T0(K)", "T0 (K)"]
    dict_measures["T_av"] = ["T_av(K)", "Taverage(K)", "T (K)"]
    dict_measures["Tp"] = ["T+(K)", "T+ (K)"]
    dict_measures["Tm"] = ["T-(K)", "T- (K)"]
    dict_measures["dTx"] = ["dTx(K)", "dTx (K)"]
    dict_measures["kxx"] = ["kxx(W/Km)", "k_xx(W/Km)", "Kxx (W / K m)"]
    dict_measures["kxy"] = ["kxy(W/mk)", "k_xy(W/Km)", "Kxy (W / K m)"]
    dict_measures["dTy"] = ["dTy(K)", "dTy (K)"]
    dict_measures["I"] = ["I(A)", "I (A)"]
    dict_measures["dTabs"] = ["dTabs", "dT_abs"]
    dict_measures["kxx/T"] = ["kxx/T"]
    dict_measures["Resistance"] = ["Resistance"]
    dict_measures["dTx/T"] = ["dTx/T"]
    dict_measures["Tp_Tm"] = ["Tp_Tm"]
    dict_measures["I_fit"] = ["I_fit"]
    dict_measures["T0_fit"] = ["T0_fit"]


    return dict_measures

def get_raw_data():
    dict_raw = dict()
    dict_raw["T0"] = ["#T0(K)"]
    dict_raw["I"] = ["I(A)"]
    dict_raw["R+_0"] = ["R+_0(V)"]
    dict_raw["R+_Q"] = ["R+_Q(V)"]
    dict_raw["R-_0"] = ["R-_0(V)"]
    dict_raw["R-_Q"] = ["R-_Q(V)"]
    dict_raw["dTy_0"] = ["dTy_0(V)"]
    dict_raw["dTy_Q"] = ["dTy_Q(V)"]
    dict_raw["dTabs_0"] = ["Tabs_0(V)"]
    dict_raw["dTabs_Q"] = ["Tabs_Q(V)"]
    dict_raw["dTx_0"] = ["dTx_0(V)"]
    dict_raw["dTx_Q"] = ["dTx_Q(V)"]

    return dict_raw

def get_parameters():
    # Creation of a dictionnary to sort other info
    dict_parameters = dict()
    dict_parameters["H"] = ["H"]
    dict_parameters["w"] = ["w"]
    dict_parameters["t"] = ["t"]
    dict_parameters["L"] = ["L"]
    dict_parameters["mount"] = ["mount"]
    dict_parameters["sample"] = ["Sample", "sample"]
    dict_parameters["date"] = ["Date", "date"]

    return dict_parameters

def get_figure_axes():
    # Creation of an internal dictionnary used to match measurements to their
    # respective axis titles to make the figures prettier.
    list_measures = list()
    list_parameters = list()
    dict_axis = dict()
    dict_axis["T_av"] = r"T ( K )"
    dict_axis["T0"] = r"$T_0$ ( K )"
    dict_axis["Tp"] = dict_axis["T_av"]
    dict_axis["Tm"] = dict_axis["T_av"]
    dict_axis["kxx"] = r"$\kappa_{\rm xx}$ ( W / K m )"
    dict_axis["dTx"] = r"$\Delta T_{\rm x}$ ( K )"
    dict_axis["dTy"] = r"$\Delta T_{\rm y}$ ( K )"
    dict_axis["kxx/T"] = r"$\kappa_{\rm xx}$/T ( W / K$^2$ m )"
    dict_axis["dTx/T"] = r"$\Delta T_{\rm x}$/T ( % )"
    dict_axis["Resistance"] = r"(T-T$_0$)/$\Delta T_{\rm x}$"
    dict_axis["kxy"] = r"$\kappa_{\rm xy}$ ( mW / K cm )"
    dict_axis["kxy/kxx"] = r"$\kappa_{\rm xy}/\kappa_{\rm xx}$ ( % )"
    dict_axis["dTy/dTx"] = r"$\Delta T_{\rm y}/\Delta T_{\rm x}$ ( % )"
    dict_axis["Tp_Tm"] = dict_axis["T_av"]
    dict_axis["T0_fit"] = dict_axis["T0"]
    dict_axis["I_fit"] = r"I ( mA )"

    return dict_axis

def get_figure_labels():
    # Same principle then before but for curve labels
    dict_labels = dict()
    dict_labels["H"] = r"H = %sT"
    dict_labels["sample"] = r"Sample: %s"
    dict_labels["date"] = r"%s"

    return dict_labels


