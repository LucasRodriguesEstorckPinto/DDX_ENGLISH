import customtkinter as ctk
import tkinter.messagebox as messagebox
import webbrowser
import sympy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re
from PIL import Image
from functools import lru_cache
from scipy.optimize import fsolve
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from sympy import Interval, Union, S, solve, log, Complement, FiniteSet, oo, Pow
from sympy.solvers.inequalities import solve_univariate_inequality
from tkinter import filedialog
from scipy.interpolate import interp1d

matplotlib.use("TkAgg")
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("green")

font = ("Segoe UI", 14)

sp.init_printing()
x = sp.symbols('x')
n = sp.symbols('n', integer=True)

idioma_atual = "PT"

dicionario_pt = {
    "Welcome to the DDX Calculator": "Bem-vindo ao DDX",
    "Open DDX Calculator": "Abrir Calculadora DDX",
    "Open DDX Manual": "Abrir Manual DDX",
    "Select Language:": "Selecione o Idioma:",
    "Domain and Range": "Domínio e Imagem",
    "Limits": "Limites",
    "Derivatives": "Derivadas",
    "Root": "Raiz",
    "Graphs": "Gráficos",
    "Implicit Derivatives": "Derivadas Implícitas",
    "Integrals": "Integral",
    "Manual": "Manual",
    "Expression:": "Expressão:",
    "Calculate": "Calcular",
    "Example": "Exemplo",
    "Function:": "Função:",
    "Point:": "Ponto:",
    "Derivative order (e.g., 1, 2, 3...):": "Ordem da derivada (ex: 1, 2, 3...):",
    "Plot Tangent (Order 1)": "Plotar Tangente (Ordem 1)",
    "Variable (empty = all):": "Variável (vazio = todas):",
    "Variable:": "Variável:",
    "Approaching:": "Tendendo a:",
    "Left": "Esquerda",
    "Right": "Direita",
    "Both": "Ambos",
    "Number:": "Número:",
    "Index:": "Índice:",
    "Function(s) (use ',' or ';'):": "Função(ões) (use ',' ou ';'):",
    "Interval(s) (use ',' or ';'):": "Intervalo(s) (use ',' ou ';'):",
    "Y Interval (optional, e.g., -5,5):": "Intervalo Y (opcional, ex: -5,5):",
    "Piecewise Function (separate with ';')": "Função por Partes (separe com ';')",
    "Show critical and inflection points": "Mostrar pontos críticos e de inflexão",
    "Plot": "Plotar",
    "Import points file": "Importar arquivo de pontos",
    "Interpolate curve": "Interpolar curva",
    "Plot imported data": "Plotar dados importados",
    "Equation (e.g., x**2 + y**2 = 1):": "Equação (ex: x**2 + y**2 = 1):",
    "Dependent Variable (e.g., y):": "Variável Dependente (ex: y):",
    "Independent Variable (e.g., x):": "Variável Independente (ex: x):",
    "Numerator Function:": "Função Numerador:",
    "Denominator Function:": "Função Denominador:",
    "Apply L'Hôpital": "Aplicar L'Hôpital",
    "Lower limit:": "Limite inferior:",
    "Upper limit:": "Limite superior:",
    "What are Domains and Ranges?": "O que são Domínios e Imagens?",
    "What is a Derivative?": "O que é uma Derivada?",
    "What are Partial Derivatives?": "O que são Derivadas Parciais?",
    "What are Limits?": "O que são Limites?",
    "When to use L'Hôpital?": "Quando usar L'Hôpital?",
    "What is an Integral?": "O que é uma Integral?",
    "Close": "Fechar"
}

def _(texto):
    if idioma_atual == "EN":
        return texto
    return dicionario_pt.get(texto, texto)

def calculo_derivadas_parciais():
    global entradafuncparcial, entradavarparcial, resultado_text_parcial
    try:
        func_str = entradafuncparcial.get()
        var_str = entradavarparcial.get().strip()
        variaveis = sorted(set(re.findall(r"[a-zA-Z]+", func_str)))
        vars_sympy = sp.symbols(" ".join(variaveis))
        expr = sp.sympify(func_str)
        resultado_text_parcial.delete("1.0", ctk.END)

        if var_str:
            var = sp.Symbol(var_str)
            derivada = sp.diff(expr, var)
            resultado_text_parcial.insert(ctk.END, f"∂f/∂{var_str} = {derivada}\n")
        else:
            for var in vars_sympy:
                derivada = sp.diff(expr, var)
                resultado_text_parcial.insert(ctk.END, f"∂f/∂{var} = {derivada}\n")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating the partial derivative.\n{e}")

def calculo_derivada():
    global resultado_text_deriv, entradaderiv, entradaponto, entradaordem
    try:
        x = sp.Symbol('x')
        func_str = entradaderiv.get()
        func = sp.sympify(func_str)
        ordem_str = entradaordem.get()
        ordem = int(ordem_str) if ordem_str else 1
        resultado_text_deriv.delete("1.0", ctk.END)
        derivada_atual = func

        for i in range(1, ordem + 1):
            derivada_atual = sp.diff(derivada_atual, x)
            resultado_text_deriv.insert(
                ctk.END,
                f"The {i}-th order derivative of the function is: {derivada_atual}\n" if idioma_atual == "EN" else f"A derivada de ordem {i} da função é: {derivada_atual}\n"
            )

        point_str = entradaponto.get()
        if point_str:
            point = float(sp.sympify(point_str))
            valor_derivada = derivada_atual.subs(x, point)
            resultado_text_deriv.insert(
                ctk.END,
                f"\nThe value of the {ordem}-th derivative at point x={point} is: {valor_derivada}\n" if idioma_atual == "EN" else f"\nO valor da derivada de ordem {ordem} no ponto x={point} é: {valor_derivada}\n"
            )
            if ordem == 1:
                coef_angular = valor_derivada
                reta = func.subs(x, point) + coef_angular * (x - point)
                resultado_text_deriv.insert(
                    ctk.END,
                    f"The equation of the tangent line is: {reta}\n\n" if idioma_atual == "EN" else f"A equação da reta tangente é: {reta}\n\n"
                )
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating the derivative: {e}")

def calculo_derivada_implicita():
    global entrada_implicita, entrada_var_dependente, entrada_var_independente, resultado_text_implicita
    try:
        eq_str = entrada_implicita.get()
        y_str = entrada_var_dependente.get().strip()
        x_str = entrada_var_independente.get().strip()

        if not eq_str or not y_str or not x_str:
            messagebox.showerror("Error", "Please fill in all fields." if idioma_atual == "EN" else "Por favor, preencha todos os campos.")
            return

        x_sym = sp.symbols(x_str)
        y_sym = sp.symbols(y_str)

        if '=' in eq_str:
            lhs_str, rhs_str = eq_str.split('=', 1)
            lhs = sp.sympify(lhs_str.strip())
            rhs = sp.sympify(rhs_str.strip())
            eq_expression = lhs - rhs
        else:
            eq_expression = sp.sympify(eq_str)

        derivada_implicita = sp.idiff(eq_expression, y_sym, x_sym)
        resultado_text_implicita.delete("1.0", ctk.END)
        resultado_text_implicita.insert(ctk.END, f"Derivative of {y_str} with respect to {x_str}:\n\n" if idioma_atual == "EN" else f"Derivada de {y_str} em relação a {x_str}:\n\n")
        resultado_text_implicita.insert(ctk.END, f"{derivada_implicita}")
    except Exception as e:
        messagebox.showerror("Error", f"Error: {e}")

def calculo_limite():
    global resultado_text_limite, entradalimit, entradavar, entradatend, direcao_var
    try:
        func_str = entradalimit.get()
        func = sp.sympify(func_str)
        variavel = sp.symbols(entradavar.get())
        valor_tendencia = float(sp.sympify(entradatend.get()))
        direcao = direcao_var.get() 

        if direcao == _("Both") or direcao == "Both":
            limite_esquerda = sp.limit(func, variavel, valor_tendencia, dir='-')
            limite_direita = sp.limit(func, variavel, valor_tendencia, dir='+')
            resultado_text_limite.delete("1.0", ctk.END)
            if limite_esquerda == limite_direita:
                resultado_text_limite.insert(ctk.END, f"The limit is: {limite_esquerda}" if idioma_atual == "EN" else f"O limite da função é: {limite_esquerda}")
            else:
                resultado_text_limite.insert(ctk.END, f"The limit does not exist." if idioma_atual == "EN" else f"O limite da função não existe.")
        else:
            dir_sym = '-' if direcao in ["Left", "Esquerda"] else '+'
            limite = sp.limit(func, variavel, valor_tendencia, dir=dir_sym)
            resultado_text_limite.delete("1.0", ctk.END)
            resultado_text_limite.insert(ctk.END, f"The limit is: {limite}" if idioma_atual == "EN" else f"O limite da função é: {limite}")
    except Exception as e:
        messagebox.showerror("Error", "Error calculating limit." if idioma_atual == "EN" else "Erro ao calcular limite.")

def raiz():
    global entradaraiz, entradaindice, resultado_text_raiz
    try:
        numero = float(entradaraiz.get())
        indice_input = entradaindice.get()
        if not indice_input:
            raise ValueError("Index not provided")
        indice = int(indice_input)

        if indice == 2:
            tolerancia = 1e-10
            x_val = numero / 2
            while True:
                raiz_value = 0.5 * (x_val + numero / x_val)
                if abs(raiz_value - x_val) < tolerancia:
                    break
                x_val = raiz_value
        else:
            raiz_value = pow(numero, 1/indice)

        resultado_text_raiz.delete("1.0", ctk.END)
        resultado_text_raiz.insert(ctk.END, f"The {indice}-th root of {numero} is: {raiz_value:.4}" if idioma_atual == "EN" else f"A raiz de índice {indice} de {numero} é: {raiz_value:.4}")
    except ValueError:
        messagebox.showerror("Error", "Invalid index." if idioma_atual == "EN" else "Índice/Número inválido.")

def validar_entrada_grafico(func_str, intervalo_str):
    if not func_str or not intervalo_str:
        raise ValueError("Empty function or interval input.")
    func_list = [f.strip() for f in func_str.split(',')]
    try:
        lower, upper = map(float, intervalo_str.split(','))
        if lower >= upper:
            raise ValueError("The lower limit must be less than the upper limit.")
    except ValueError as e:
        raise ValueError(f"Invalid interval: {str(e)}")
    for f in func_list:
        try:
            sp.sympify(f)
        except sp.SympifyError:
            raise ValueError(f"Invalid function: {f}")
    return func_list, lower, upper

@lru_cache(maxsize=128)
def calcular_derivadas(func, x):
    fprime = sp.diff(func, x)
    fsecond = sp.diff(fprime, x)
    return fprime, fsecond

def encontrar_assintota_obliqua(func, x):
    numer, denom = func.as_numer_denom()
    deg_numer = sp.degree(numer, gen=x)
    deg_denom = sp.degree(denom, gen=x)
    if deg_numer - deg_denom == 1:
        coef = sp.limit(func/denom, x, sp.oo)
        intercept = sp.limit(func - coef*x, x, sp.oo)
        return coef, intercept
    return None, None

def numerical_roots(expr, var, a, b, num_points=500):
    x_vals = np.linspace(a, b, num_points)
    roots = []
    try:
        expr_func = sp.lambdify(var, expr, 'numpy')
        y_vals = expr_func(x_vals)
        if np.isscalar(y_vals):
            y_vals = np.full_like(x_vals, y_vals)
        for i in range(len(x_vals)-1):
            try:
                val1, val2 = float(y_vals[i]), float(y_vals[i+1])
                if np.isfinite(val1) and np.isfinite(val2):
                    if np.sign(val1) * np.sign(val2) < 0:
                        def eq_func(v):
                            res = expr_func(v[0]) 
                            return float(res) if np.isfinite(float(res)) else 0.0
                        root_array = fsolve(eq_func, (x_vals[i] + x_vals[i+1]) / 2)
                        root = float(root_array[0])
                        if a <= root <= b and not any(abs(root - r) < 1e-6 for r in roots):
                            roots.append(root)
                    elif val1 == 0.0:
                        root = float(x_vals[i])
                        if a <= root <= b and not any(abs(root - r) < 1e-6 for r in roots):
                            roots.append(root)
            except Exception:
                continue 
    except Exception as e:
        print(f"Erro no rastreio numérico de raízes: {e}")
    return sorted(roots)

def ajustar_amostragem(lower, upper, num_points_base=200):
    if upper - lower > 100:
        return np.linspace(lower, upper, num_points_base)
    return np.linspace(lower, upper, min(num_points_base * 2, 800))

def carregar_arquivo_pontos():
    global dados_x, dados_y, botao_plot_dados, check_interpolar, resultado_text_grafico
    file_path = filedialog.askopenfilename(filetypes=[("Text Files", "*.txt")])
    if not file_path:
        return
    try:
        with open(file_path, 'r') as file:
            linhas = file.readlines()
            dados_lidos = []
            for linha in linhas:
                linha = re.sub(r'[,\t]+', ' ', linha.strip()) 
                partes = linha.split()
                if len(partes) >= 2:
                    dados_lidos.append([float(partes[0]), float(partes[1])])
        dados = np.array(dados_lidos)
        dados_x = dados[:, 0]
        dados_y = dados[:, 1]
        if len(dados_x) < 2:
            raise ValueError("Minimum of two points required.")
        check_interpolar.pack(pady=5, anchor="w")
        botao_plot_dados.pack(pady=5, padx=5, anchor="w")
        resultado_text_grafico.insert(ctk.END, "\nFile loaded successfully.\n" if idioma_atual == "EN" else "\nArquivo carregado com sucesso.\n")
    except Exception as e:
        messagebox.showerror("Error", f"Error: {str(e)}")

def plotar_dados_importados():
    global dados_x, dados_y, grafico_canvas, grafico_toolbar, frame_grafico_container
    if dados_x is None or dados_y is None:
        return
    try:
        if frame_grafico_container:
            for widget in frame_grafico_container.winfo_children():
                widget.destroy()
        ordenados = sorted(zip(dados_x, dados_y), key=lambda par: par[0])
        x_ord, y_ord = zip(*ordenados)
        x_ord = np.array(x_ord)
        y_ord = np.array(y_ord)
        fig, ax = plt.subplots(figsize=(10, 6))

        if interpolar_var.get():
            if len(x_ord) >= 4:
                try:
                    f_interp = interp1d(x_ord, y_ord, kind='cubic')
                    x_interp = np.linspace(x_ord[0], x_ord[-1], 500)
                    y_interp = f_interp(x_interp)
                    ax.plot(x_interp, y_interp, label="Interpolated", linewidth=2.5, color='cyan', zorder=4)
                except Exception as e:
                    pass
        ax.scatter(x_ord, y_ord, color='red', s=60, zorder=5, label="Points")
        ax.axhline(0, color='black', lw=1.2, linestyle='dashed')
        ax.axvline(0, color='black', lw=1.2, linestyle='dashed')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title("Data plot")
        ax.legend()
        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frame_grafico_container)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        grafico_canvas = canvas
        toolbar = NavigationToolbar2Tk(canvas, frame_grafico_container)
        toolbar.update()
        toolbar.pack()
        grafico_toolbar = toolbar
    except Exception as e:
        messagebox.showerror("Error", f"Error: {str(e)}")

def plot_grafico():
    global resultado_text_grafico, entrada_grafico, intervalo, entrada_intervalo_y, show_points_var, is_piecewise_var, grafico_canvas, grafico_toolbar, frame_grafico_container
    try:
        if frame_grafico_container:
            for widget in frame_grafico_container.winfo_children():
                widget.destroy()
        plt.style.use('ggplot')
        plt.rcParams.update({
            'font.size': 12,
            'axes.titlesize': 18,
            'axes.labelsize': 14,
            'legend.fontsize': 12
        })
        fig, ax = plt.subplots(figsize=(10, 6))
        result_text = ""
        func_input_str = entrada_grafico.get()
        intervalo_input_str = intervalo.get()
        lista_de_pedacos = []
        full_lower, full_upper = float('inf'), float('-inf')

        if is_piecewise_var.get():
            func_list_str = [f.strip() for f in func_input_str.split(';') if f.strip()]
            interval_list_str = [i.strip() for i in intervalo_input_str.split(';') if i.strip()]
            if not func_list_str or not interval_list_str:
                raise ValueError("Insert functions and intervals.")
            if len(func_list_str) != len(interval_list_str):
                raise ValueError("Mismatch in number of functions and intervals.")
            for f_str, i_str in zip(func_list_str, interval_list_str):
                func_sym = sp.sympify(f_str)
                parts = [p.strip() for p in i_str.split(',')]
                lower = float(sp.N(sp.sympify(parts[0])))
                upper = float(sp.N(sp.sympify(parts[1])))
                if lower > upper:
                    lower, upper = upper, lower
                lista_de_pedacos.append( (func_sym, lower, upper) )
                full_lower = min(full_lower, lower)
                full_upper = max(full_upper, upper)
        else:
            func_list_str = [f.strip() for f in func_input_str.split(',') if f.strip()]
            parts = [p.strip() for p in intervalo_input_str.split(',')]
            lower = float(sp.N(sp.sympify(parts[0])))
            upper = float(sp.N(sp.sympify(parts[1])))
            if lower > upper:
                lower, upper = upper, lower
            full_lower, full_upper = lower, upper
            for f_str in func_list_str:
                func_sym = sp.sympify(f_str)
                lista_de_pedacos.append( (func_sym, lower, upper) )

        y_lower = y_upper = None
        try:
            y_intervalo_str = entrada_intervalo_y.get().strip()
            if y_intervalo_str:
                parts = [p.strip() for p in y_intervalo_str.split(',')]
                if len(parts) == 2:
                    y_lower = float(sp.N(sp.sympify(parts[0])))
                    y_upper = float(sp.N(sp.sympify(parts[1])))
                    if y_lower > y_upper:
                        y_lower, y_upper = y_upper, y_lower
        except Exception:
            y_lower = y_upper = None

        for i, (func_sym, lower, upper) in enumerate(lista_de_pedacos):
            func_numeric = sp.lambdify(x, func_sym, 'numpy')
            x_vals = ajustar_amostragem(lower, upper)
            y_vals = func_numeric(x_vals)
            ax.plot(x_vals, y_vals, label=f'${sp.latex(func_sym)}$ in [{lower},{upper}]', linewidth=2.5, color=f'C{i}')

        ax.axhline(0, color='black', lw=1.2, linestyle='dashed', zorder=3)
        ax.axvline(0, color='black', lw=1.2, linestyle='dashed', zorder=3)
        ax.set_xlabel('x', fontsize=14)
        ax.set_ylabel('y', fontsize=14)
        ax.set_xlim(full_lower, full_upper)
        if y_lower is not None and y_upper is not None:
            ax.set_ylim(y_lower, y_upper)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frame_grafico_container)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        grafico_canvas = canvas
        toolbar = NavigationToolbar2Tk(canvas, frame_grafico_container)
        toolbar.update()
        toolbar.pack()
        grafico_toolbar = toolbar
    except Exception as e:
        messagebox.showerror("Error", f"Error: {str(e)}")

def validar_entrada(func_str):
    pattern = r'^[a-zA-Z0-9\s+\-*/().^sincoslogexp]+$'
    if not re.match(pattern, func_str):
        raise ValueError("Invalid input.")
    return func_str.replace("^", "**").replace("sen", "sin").replace("arctg", "atan").replace("arcsen", "asin").replace("arccos", "acos")

def formatar_conjunto(conjunto):
    if isinstance(conjunto, str):
        return conjunto
    try:
        return str(sp.simplify(conjunto))
    except Exception:
        return repr(conjunto)

def calcular_dominio(func, x):
    try:
        func = sp.sympify(func)
        dominio = S.Reals
        return dominio
    except Exception as e:
        return f"Error: {e}"

def calculo_dominio_imagem():
    global resultado_text_dom, entradadom, grafico_label, x
    try:
        func_str = entradadom.get()
        func_str = validar_entrada(func_str) 
        func = sp.sympify(func_str)
        dominio = calcular_dominio(func, x)
        resultado = f"Domain: {formatar_conjunto(dominio)}"
        resultado_text_dom.delete("1.0", ctk.END)
        resultado_text_dom.insert(ctk.END, resultado)
    except Exception as e:
        messagebox.showerror("Error", f"Error: {str(e)}")

def calculo_integral():
    global resultado_text_integral, entrada_integrais, entrada_limite_inf, entrada_limite_sup
    try:
        func_str = entrada_integrais.get()
        x = sp.symbols('x')
        func = sp.sympify(func_str)
        limite_inf_str = entrada_limite_inf.get().strip()
        limite_sup_str = entrada_limite_sup.get().strip()
        if limite_inf_str and limite_sup_str:
            limite_inf = float(sp.sympify(limite_inf_str))
            limite_sup = float(sp.sympify(limite_sup_str))
            integral_def = sp.integrate(func, (x, limite_inf, limite_sup))
            resultado_text_integral.delete("1.0", ctk.END)
            resultado_text_integral.insert(ctk.END, f"Integral: {integral_def}\n")
        else:
            integral = sp.integrate(func, x)
            resultado_text_integral.delete("1.0", ctk.END)
            resultado_text_integral.insert(ctk.END, f"Integral: {integral} + C\n")
    except Exception as e:
        messagebox.showerror("Error", f"Error: {e}")

def plot_func_tangente():
    try:
        x = sp.Symbol('x')
        func_str = entradaderiv.get()
        func = sp.sympify(func_str)
        point = float(sp.sympify(entradaponto.get()))
        derivada = sp.diff(func, x)
        coef_angular = derivada.subs(x, point)
        reta = func.subs(x, point) + coef_angular * (x - point)
        func_num = sp.lambdify(x, func, "numpy")
        reta_num = sp.lambdify(x, reta, "numpy")
        x_vals = np.linspace(point - 10, point + 10, 400)
        plt.figure()
        y_func = func_num(x_vals)
        if np.isscalar(y_func):
            y_func = np.full_like(x_vals, y_func)
        y_reta = reta_num(x_vals)
        if np.isscalar(y_reta):
            y_reta = np.full_like(x_vals, y_reta)
        plt.plot(x_vals, y_func, label=f"f(x) = {func_str}")
        plt.plot(x_vals, y_reta, label=f"Tangent at x = {point}")
        plt.axhline(0, color='red', lw=0.8)
        plt.axvline(0, color='red', lw=0.8)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.grid(True)
        plt.show()
    except Exception as e:
        messagebox.showerror("Error", f"Error: {e}")

def aplicar_lhopital(f_str, g_str, ponto_str, direcao='Both'):
    from sympy import limit, sympify, diff, simplify
    f = sympify(f_str)
    g = sympify(g_str)
    ponto = sympify(ponto_str)
    passos = []
    def calcular(expr, lado):
        return limit(expr, x, ponto, dir=lado)
    lados = ["+"] if direcao in ["+"] else ["-"] if direcao == "-" else ["+", "-"]
    for lado in lados:
        passos.append(f"Analyzing lateral limit: {'right' if lado == '+' else 'left'}")
        try:
            lim_f = calcular(f, lado)
            lim_g = calcular(g, lado)
            passos.append(f"  lim(x→{ponto}{lado}) {f} = {lim_f}")
            passos.append(f"  lim(x→{ponto}{lado}) {g} = {lim_g}")
        except Exception as e:
            passos.append(f"  Error calculating limits: {e}")
            continue
        formas_validas = [
            abs(lim_f.evalf()) < 1e-10 and abs(lim_g.evalf()) < 1e-10,
            lim_f.is_infinite and lim_g.is_infinite]
        if not any(formas_validas):
            passos.append("  ❌ L'Hôpital's Rule DOES NOT apply — form is not indeterminate.")
            continue
        passos.append("  ✅ Indeterminate form detected. Applying L'Hôpital:")
        i = 1
        num, den = f, g
        while i <= 10:
            num_deriv = diff(num, x)
            den_deriv = diff(den, x)
            passos.append(f"    Iteration {i}:")
            passos.append(f"      f'(x) = {num_deriv}")
            passos.append(f"      g'(x) = {den_deriv}")
            try:
                lim_num = calcular(num_deriv, lado)
                lim_den = calcular(den_deriv, lado)
                passos.append(f"      lim(x→{ponto}{lado}) f'(x) = {lim_num}")
                passos.append(f"      lim(x→{ponto}{lado}) g'(x) = {lim_den}")
            except Exception as e:
                passos.append(f"      Error calculating limits of derivatives: {e}")
                break
            if lim_den != 0 and lim_den.is_number and lim_num.is_number:
                resultado = simplify(lim_num / lim_den)
                passos.append(f"      ✅ Final result after {i} iteration(s): {resultado}")
                break
            num, den = num_deriv, den_deriv
            i += 1
            if i > 10:
                passos.append("      ❌ Maximum number of iterations reached.")
                break
        passos.append("")
    return passos

def calculo_lhopital():
    global entrada_num, entrada_den, entrada_ponto, direcao_lhopital, resultado_text_lhopital
    try:
        num = entrada_num.get()
        den = entrada_den.get()
        ponto = entrada_ponto.get()
        direcao = direcao_lhopital.get()
        passos = aplicar_lhopital(num, den, ponto, direcao)
        resultado_text_lhopital.delete("1.0", ctk.END)
        for passo in passos:
            resultado_text_lhopital.insert(ctk.END, passo + "\n")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def exemplo_lhopital():
    entrada_num.delete(0, ctk.END)
    entrada_den.delete(0, ctk.END)
    entrada_ponto.delete(0, ctk.END)
    entrada_num.insert(0, "sin(x)")
    entrada_den.insert(0, "x")
    entrada_ponto.insert(0, "0")

def exemplo_raiz():
    if idioma_atual == "PT":
        example_text = ("Exemplo de Raiz Quadrada:\n"
            "Número: 256\n"
            "Definição: A raiz quadrada de um número é um valor que, quando multiplicado por si mesmo, "
            "resulta no número original.\n"
            "Cálculo: A raiz quadrada de 256 é 16, porque 16 * 16 = 256.\n"
            "Propriedades: A raiz quadrada de um número positivo é sempre um número positivo. "
            "Neste caso, a raiz quadrada de 256 é um valor inteiro exato, 16.")
    else:
        example_text = ("Example of Square Root:\n"
            "Number: 256\n"
            "Definition: The square root of a number is a value that, when multiplied by itself, "
            "results in the original number.\n"
            "Calculation: The square root of 256 is 16, because 16 * 16 = 256.\n"
            "Properties: The square root of a positive number is always a positive number. "
            "In this case, the square root of 256 is an exact integer value, 16.")
    resultado_text_raiz.delete("1.0", ctk.END)
    resultado_text_raiz.insert(ctk.END, example_text)

def exemplo_dominio_imagem():
    if idioma_atual == "PT":
        example_text = (
            "Exemplo de Domínio e Imagem:\n"
            "Função: f(x) = 1/(x-2)\n"
            "Domínio: Todos os valores de x, exceto x=2. Isso ocorre porque a função se torna indefinida quando x=2, "
            "pois resultaria em uma divisão por zero.\n"
            "Imagem: Todos os valores reais, exceto f(x)=0. A função nunca toca o eixo x, "
            "pois não há valor de x que torne a função igual a zero."
        )
    else:
        example_text = (
            "Example of Domain and Range:\n"
            "Function: f(x) = 1/(x-2)\n"
            "Domain: All values of x, except x=2. This is because the function becomes undefined when x=2, "
            "as it would result in a division by zero.\n"
            "Range: All real values, except f(x)=0. The function never touches the x-axis, "
            "as there is no value of x that makes the function equal to zero."
        )
    resultado_text_dom.delete("1.0", ctk.END)
    resultado_text_dom.insert(ctk.END, example_text)

def exemplo_limite():
    if idioma_atual == "PT":
        example_text = (
            "Exemplo de Limite:\n"
            "Função: f(x) = (x^2 - 1)/(x - 1)\n"
            "Para calcular o limite de f(x) quando x se aproxima de 1, simplificamos a função:\n"
            "f(x) = (x + 1) para x ≠ 1.\n"
            "Então, o limite de f(x) quando x se aproxima de 1 é 2.\n"
            "Lembre-se que o limite se refere ao valor do qual a função se aproxima."
        )
    else:
        example_text = (
            "Example of Limit:\n"
            "Function: f(x) = (x^2 - 1)/(x - 1)\n"
            "To calculate the limit of f(x) as x approaches 1, we simplify the function:\n"
            "f(x) = (x + 1) for x ≠ 1.\n"
            "Then, the limit of f(x) as x approaches 1 is 2.\n"
            "Remember that the limit refers to the value the function approaches as x approaches 1."
        )
    resultado_text_limite.delete("1.0", ctk.END)
    resultado_text_limite.insert(ctk.END, example_text)

def exemplo_derivada():
    if idioma_atual == "PT":
        example_text = (
            "Exemplo de Derivada e Tangente:\n"
            "Função: f(x) = x^2\n"
            "Derivada: f'(x) = 2x. Isso representa a inclinação da função em qualquer ponto x.\n"
            "No ponto x=3, f'(3) = 6. Isso significa que a inclinação da tangente à curva no ponto (3, f(3)) é 6.\n"
            "A equação da reta tangente é dada por: y = f(3) + f'(3)*(x - 3)\n"
            "Neste caso, a reta tangente é y = 9 + 6(x - 3), simplificando: y = 6x - 9."
        )
    else:
        example_text = (
            "Example of Derivative and Tangent:\n"
            "Function: f(x) = x^2\n"
            "Derivative: f'(x) = 2x. This represents the slope of the function at any point x.\n"
            "At point x=3, f'(3) = 6. This means the slope of the tangent to the curve at point (3, f(3)) is 6.\n"
            "The equation of the tangent line is given by: y = f(3) + f'(3)*(x - 3)\n"
            "In this case, the tangent line is y = 9 + 6(x - 3), simplifying: y = 6x - 9."
        )
    resultado_text_deriv.delete("1.0", ctk.END)
    resultado_text_deriv.insert(ctk.END, example_text)

def exemplo_derivada_implicita():
    entrada_implicita.delete(0, ctk.END)
    entrada_var_dependente.delete(0, ctk.END)
    entrada_var_independente.delete(0, ctk.END)
    entrada_implicita.insert(0, "x**2 + y**2 = 1")
    entrada_var_dependente.insert(0, "y")
    entrada_var_independente.insert(0, "x")

def exemplo_derivada_parcial():
    entradafuncparcial.delete(0, ctk.END)
    entradavarparcial.delete(0, ctk.END)
    entradafuncparcial.insert(0, "x**2 * y + sin(z)")
    entradavarparcial.insert(0, "x")

def exemplo_integral():
    if idioma_atual == "PT":
        example_text = (
            "Exemplo de Integral:\n"
            "Função: f(x) = x^2\n"
            "Integral Indefinida: ∫x^2 dx = (1/3)x^3 + C, onde C é a constante de integração.\n"
            "Integral Definida de 0 a 2: ∫(de 0 a 2) x^2 dx = [(1/3)x^3] de 0 a 2 = (8/3) - 0 = 8/3.\n"
            "Isso representa a área sob a curva de f(x) entre x=0 e x=2."
        )
    else:
        example_text = (
            "Example of Integral:\n"
            "Function: f(x) = x^2\n"
            "Indefinite Integral: ∫x^2 dx = (1/3)x^3 + C, where C is the constant of integration.\n"
            "Definite Integral from 0 to 2: ∫(from 0 to 2) x^2 dx = [(1/3)x^3] from 0 to 2 = (8/3) - 0 = 8/3.\n"
            "This represents the area under the curve of f(x) between x=0 and x=2."
        )
    resultado_text_integral.delete("1.0", ctk.END)
    resultado_text_integral.insert(ctk.END, example_text)

# ======== FUNÇÕES DE EXPLICAÇÃO RESTAURADAS ========

def abrir_explicacao_integral():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title(_("What is an Integral?"))
    janela_explicacao.geometry("500x300")
    janela_explicacao.lift()

    if idioma_atual == "PT":
        texto_explicacao = "A integral de uma função representa a área sob a curva dessa função em um dado intervalo.\nÉ usada para calcular áreas, volumes e resolver problemas físicos, como trabalho e deslocamento.\n\nFonte: Stewart, James. Cálculo. 8ª edição."
    else:
        texto_explicacao = "The integral of a function represents the area under the curve of that function in a given interval.\nIt is used to calculate areas, volumes, and solve physical problems such as work and displacement.\n\nSource: Stewart, James. Calculus. 8th edition."

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text=_("Close"), command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_derivada():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title(_("What is a Derivative?"))
    janela_explicacao.geometry("500x300")
    janela_explicacao.lift()

    if idioma_atual == "PT":
        texto_explicacao = "A derivada de uma função representa a taxa de variação dessa função em um determinado ponto.\nÉ usada para calcular velocidades, acelerações e resolver problemas físicos, como otimização e crescimento populacional.\n\nFonte: Stewart, James. Cálculo. 8ª edição."
    else:
        texto_explicacao = "The derivative of a function represents the rate of change of that function at a given point.\nIt is used to calculate velocities, accelerations, and solve physical problems such as optimization and population growth.\n\nSource: Stewart, James. Calculus. 8th edition."

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text=_("Close"), command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_limites():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title(_("What are Limits?"))
    janela_explicacao.geometry("500x300")
    janela_explicacao.lift()

    if idioma_atual == "PT":
        texto_explicacao = "O limite de uma função descreve o comportamento dessa função quando a variável independente se aproxima de um determinado valor.\nÉ usado para definir derivadas, integrais e resolver problemas envolvendo continuidade e comportamento assintótico.\n\nFonte: Stewart, James. Cálculo. 8ª edição."
    else:
        texto_explicacao = "The limit of a function describes the behavior of that function as the independent variable approaches a certain value.\nIt is used to define derivatives, integrals, and solve problems involving continuity and asymptotic behavior.\n\nSource: Stewart, James. Calculus. 8th edition."

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text=_("Close"), command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_dominios():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title(_("What are Domains and Ranges?"))
    janela_explicacao.geometry("500x320")
    janela_explicacao.lift()

    if idioma_atual == "PT":
        texto_explicacao = "O domínio de uma função é o conjunto de todos os valores de entrada para os quais a função está definida.\nA imagem de uma função é o conjunto de todos os valores de saída que a função pode assumir.\nSão usados para entender o comportamento e as restrições de funções em diversos contextos matemáticos e aplicados.\n\nFonte: Stewart, James. Cálculo. 8ª edição."
    else:
        texto_explicacao = "The domain of a function is the set of all input values for which the function is defined.\nThe range of a function is the set of all output values that the function can take.\nThey are used to understand the behavior and restrictions of functions in various mathematical and applied contexts.\n\nSource: Stewart, James. Calculus. 8th edition."

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text=_("Close"), command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_lhopital():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title(_("When to use L'Hôpital?"))
    janela_explicacao.geometry("600x350")
    janela_explicacao.lift()

    if idioma_atual == "PT":
        texto_explicacao = "A Regra de L'Hôpital é usada para resolver limites que apresentam formas indeterminadas, como 0/0 ou ∞/∞.\n\nSejam f(x) e g(x) funções diferenciáveis num intervalo aberto contendo 'a', e se:\n- lim(x→a) f(x) = 0 e lim(x→a) g(x) = 0  ou\n- lim(x→a) f(x) = ∞ e lim(x→a) g(x) = ∞\n\nEntão:\n    lim(x→a) f(x)/g(x) = lim(x→a) f'(x)/g'(x)\ndesde que o limite da derivada exista.\n\nA regra pode ser aplicada repetidamente até que a indeterminação desapareça."
    else:
        texto_explicacao = "L'Hôpital's Rule is used to solve limits that present indeterminate forms,\nsuch as 0/0 or ∞/∞.\n\nLet f(x) and g(x) be differentiable functions on an open interval containing 'a', and if:\n- lim(x→a) f(x) = 0 and lim(x→a) g(x) = 0  or\n- lim(x→a) f(x) = ∞ and lim(x→a) g(x) = ∞\n\nThen:\n    lim(x→a) f(x)/g(x) = lim(x→a) f'(x)/g'(x)\nprovided this derivative limit exists.\n\nThe rule can be applied repeatedly until the indeterminacy disappears."

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, justify="left", wraplength=580, font=("Segoe UI", 14))
    label_texto.pack(padx=20, pady=20)
    
    botao_fechar = ctk.CTkButton(janela_explicacao, text=_("Close"), command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_derivadas_parciais():
    pass

class ModernEntry(ctk.CTkEntry):
    def get(self):
        text = super().get()
        text = re.sub(r'(?<=\d)(?=pi\b)', '*', text, flags=re.IGNORECASE)
        text = re.sub(r'(?<=\d)(?=e\b)', '*', text, flags=re.IGNORECASE)
        text = re.sub(r'\bpi\b', 'pi', text, flags=re.IGNORECASE)
        text = re.sub(r'\be\b', 'E', text, flags=re.IGNORECASE)
        return text

def labeled_input(parent, label_text):
    frame = ctk.CTkFrame(parent)
    frame.pack(anchor="w", pady=5, padx=5, fill="x")
    label = ctk.CTkLabel(frame, text=_(label_text), font=font)
    label.pack(anchor="w", padx=5)
    entry = ctk.CTkEntry(frame, width=400, height=30, corner_radius=5, font=font)
    entry.pack(padx=5, pady=5, anchor="w")
    return entry

def botao(parent, func, texto):
    btn = ctk.CTkButton(parent, text=_(texto), command=func, width=200, corner_radius=5, font=font)
    btn.pack(pady=5, padx=5, anchor="w")

def validar_expressao_em_tempo_real(entry_widget):
    try:
        expr = entry_widget.get()
        sp.sympify(expr)
        entry_widget.configure(border_color="green")
    except Exception:
        entry_widget.configure(border_color="red")

def aplicar_validacao_em_tempo_real(entry_widget):
    var = ctk.StringVar()
    entry_widget.configure(textvariable=var)
    var.trace_add("write", lambda *args: validar_expressao_em_tempo_real(entry_widget))


# ====================== TELA INICIAL =========================
class InitialPage(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Home Page")
        self.geometry("500x400")
        self.resizable(False, False)
        self.configure(padx=20, pady=20)

        self.label_titulo = ctk.CTkLabel(self, text=_("Welcome to the DDX Calculator"), font=("Segoe UI", 20, "bold"))
        self.label_titulo.pack(pady=20)

        ctk.CTkLabel(self, text="Select Language / Selecione o Idioma:", font=("Segoe UI", 14)).pack(pady=5)
        
        self.lang_var = ctk.StringVar(value=idioma_atual)
        self.lang_menu = ctk.CTkOptionMenu(
            self, 
            variable=self.lang_var, 
            values=["EN", "PT"], 
            command=self.mudar_idioma
        )
        self.lang_menu.pack(pady=10)

        self.open_calculator_btn = ctk.CTkButton(self, text=_("Open DDX Calculator"), command=self.open_calculator, width=250)
        self.open_calculator_btn.pack(pady=10)

        self.manual_btn = ctk.CTkButton(
            self,
            text=_("Open DDX Manual"),
            command=lambda: webbrowser.open('https://docs.google.com/document/d/1hvcUL36juGBm_8lsdOpPrMLWzmYnGvakKHaMj1BbxlY/edit?usp=sharing'),
            width=250
        )
        self.manual_btn.pack(pady=10)

    def mudar_idioma(self, novo_idioma):
        global idioma_atual
        idioma_atual = novo_idioma
        self.label_titulo.configure(text=_("Welcome to the DDX Calculator"))
        self.open_calculator_btn.configure(text=_("Open DDX Calculator"))
        self.manual_btn.configure(text=_("Open DDX Manual"))

    def open_calculator(self):
        self.destroy()
        app = App()
        app.mainloop()


# ====================== APLICAÇÃO PRINCIPAL =========================
class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("DDX Calculator")
        self.geometry("1400x800")
        self.minsize(1000, 700)

        self.create_widgets()

    def create_widgets(self):
        tabview = ctk.CTkTabview(self)
        tabview.pack(padx=10, pady=10, fill="both", expand=True)

        abas_originais = ["Domain and Range", "Limits", "Derivatives", "Root", "Graphs", "Implicit Derivatives", "L'Hôpital", "Integrals", "Manual"]
        frames = {}
        for aba in abas_originais:
            nome_traduzido = _(aba)
            frames[aba] = tabview.add(nome_traduzido)

        self.aba_dominio(frames["Domain and Range"])
        self.aba_limites(frames["Limits"])
        self.aba_derivadas(frames["Derivatives"])
        self.aba_raiz(frames["Root"])
        self.aba_graficos(frames["Graphs"])
        self.aba_derivadas_implicitas(frames["Implicit Derivatives"])
        self.aba_lhopital(frames["L'Hôpital"])
        self.aba_integrais(frames["Integrals"])
        self.aba_manual(frames["Manual"])

    # ====================== ESTRUTURA PADRÃO DAS ABAS =========================
    def estrutura_aba(self, frame):
        container = ctk.CTkFrame(frame)
        container.pack(fill="both", expand=True, padx=20, pady=20)

        left_frame = ctk.CTkFrame(container)
        left_frame.pack(side="left", fill="y", padx=10, pady=10)

        right_frame = ctk.CTkFrame(container)
        right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

        return left_frame, right_frame

    # ====================== ABA DOMÍNIO =========================
    def aba_dominio(self, frame):
        global entradadom, resultado_text_dom
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("What are Domains and Ranges?"), command=abrir_explicacao_dominios).pack(pady=5, anchor="w")

        entradadom = labeled_input(left, "Expression:")
        aplicar_validacao_em_tempo_real(entradadom)
        botao(left, calculo_dominio_imagem, "Calculate")
        botao(left, exemplo_dominio_imagem, "Example")

        resultado_text_dom = ctk.CTkTextbox(right, font=font)
        resultado_text_dom.pack(fill="both", expand=True)

    # ====================== ABA DERIVADAS =========================
    def aba_derivadas(self, frame):
        global entradaderiv, entradaponto, entradaordem, resultado_text_deriv
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("What is a Derivative?"), command=abrir_explicacao_derivada).pack(pady=5, anchor="w")

        entradaderiv = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entradaderiv)
        entradaponto = labeled_input(left, "Point:")
        entradaordem = labeled_input(left, "Derivative order (e.g., 1, 2, 3...):")

        botao(left, calculo_derivada, "Calculate")
        botao(left, exemplo_derivada, "Example")
        botao(left, plot_func_tangente, "Plot Tangent (Order 1)")

        resultado_text_deriv = ctk.CTkTextbox(right, font=font)
        resultado_text_deriv.pack(fill="both", expand=True)


        try:
            img = ctk.CTkImage(Image.open("deriva.png"), size=(250, 120))
            ctk.CTkLabel(right, image=img, text="").pack(pady=10)
        except Exception:
            pass

    # ====================== ABA DERIVADAS PARCIAIS =========================
    def aba_derivadas_parciais(self, frame):
        global entradafuncparcial, entradavarparcial, resultado_text_parcial
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("What are Partial Derivatives?"), command=abrir_explicacao_derivadas_parciais).pack(pady=5, anchor="w")

        entradafuncparcial = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entradafuncparcial)
        entradavarparcial = labeled_input(left, "Variable (empty = all):")

        botao(left, calculo_derivadas_parciais, "Calculate")
        botao(left, exemplo_derivada_parcial, "Example")

        resultado_text_parcial = ctk.CTkTextbox(right, font=font)
        resultado_text_parcial.pack(fill="both", expand=True)


    # ====================== ABA LIMITES =========================
    def aba_limites(self, frame):
        global entradalimit, entradavar, entradatend, direcao_var, resultado_text_limite
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("What are Limits?"), command=abrir_explicacao_limites).pack(pady=5, anchor="w")

        entradalimit = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entradalimit)
        entradavar = labeled_input(left, "Variable:")
        entradatend = labeled_input(left, "Approaching:")
        aplicar_validacao_em_tempo_real(entradatend)

        direcao_var = ctk.StringVar(value=_("Both"))
        ctk.CTkOptionMenu(left, variable=direcao_var, values=[_("Left"), _("Right"), _("Both")]).pack(pady=5, anchor="w")

        botao(left, calculo_limite, "Calculate")
        botao(left, exemplo_limite, "Example")

        resultado_text_limite = ctk.CTkTextbox(right, font=font)
        resultado_text_limite.pack(fill="both", expand=True)

        try:
            img = ctk.CTkImage(Image.open("limit.png"), size=(250, 120))
            ctk.CTkLabel(right, image=img, text="").pack(pady=10)
        except Exception:
            pass

    # ====================== ABA RAIZ =========================
    def aba_raiz(self, frame):
        global entradaraiz, entradaindice, resultado_text_raiz
        left, right = self.estrutura_aba(frame)

        entradaraiz = labeled_input(left, "Number:")
        aplicar_validacao_em_tempo_real(entradaraiz)
        entradaindice = labeled_input(left, "Index:")
        aplicar_validacao_em_tempo_real(entradaindice)

        botao(left, raiz, "Calculate")
        botao(left, exemplo_raiz, "Example")

        resultado_text_raiz = ctk.CTkTextbox(right, font=font)
        resultado_text_raiz.pack(fill="both", expand=True)

      
        try:
            nome_imagem = "raiz2.png" if idioma_atual == "PT" else "raiz.png"
            img = ctk.CTkImage(Image.open(nome_imagem), size=(250, 120))
            ctk.CTkLabel(right, image=img, text="").pack(pady=10)
        except Exception:
            pass

    # ====================== ABA GRÁFICOS =========================
    def aba_graficos(self, frame):
        global entrada_grafico, intervalo, entrada_intervalo_y, show_points_var, resultado_text_grafico, frame_grafico_container, is_piecewise_var
        global interpolar_var, botao_plot_dados, check_interpolar, font 

        left, right = self.estrutura_aba(frame)

        entrada_grafico = labeled_input(left, "Function(s) (use ',' or ';'):")
        aplicar_validacao_em_tempo_real(entrada_grafico)

        intervalo = labeled_input(left, "Interval(s) (use ',' or ';'):")
        aplicar_validacao_em_tempo_real(intervalo)

        entrada_intervalo_y = labeled_input(left, "Y Interval (optional, e.g., -5,5):")
        aplicar_validacao_em_tempo_real(entrada_intervalo_y)

        is_piecewise_var = ctk.BooleanVar(value=False, master=left)
        ctk.CTkCheckBox(left, text=_("Piecewise Function (separate with ';')"), variable=is_piecewise_var).pack(pady=5, anchor="w")

        show_points_var = ctk.BooleanVar(value=False, master=left)
        ctk.CTkCheckBox(left, text=_("Show critical and inflection points"), variable=show_points_var).pack(pady=5, anchor="w")

        botao(left, plot_grafico, "Plot")

        ctk.CTkButton(left, text=_("Import points file"), command=carregar_arquivo_pontos).pack(pady=10, anchor="w")
        interpolar_var = ctk.BooleanVar(value=False,master=left)
        check_interpolar = ctk.CTkCheckBox(left, text=_("Interpolate curve"), variable=interpolar_var)
        botao_plot_dados = ctk.CTkButton(left, text=_("Plot imported data"), command=plotar_dados_importados)

        resultado_text_grafico = ctk.CTkTextbox(right, font=font, height=150)
        resultado_text_grafico.pack(fill="x", pady=(0, 10))

        frame_grafico_container = ctk.CTkFrame(right)
        frame_grafico_container.pack(fill="both", expand=True)

    # ====================== ABA DERIVADAS IMPLÍCITAS =========================
    def aba_derivadas_implicitas(self, frame):
        global entrada_implicita, entrada_var_dependente, entrada_var_independente, resultado_text_implicita
        left, right = self.estrutura_aba(frame)

        entrada_implicita = labeled_input(left, "Equation (e.g., x**2 + y**2 = 1):")
        aplicar_validacao_em_tempo_real(entrada_implicita)
        entrada_var_dependente = labeled_input(left, "Dependent Variable (e.g., y):")
        entrada_var_independente = labeled_input(left, "Independent Variable (e.g., x):")

        botao(left, calculo_derivada_implicita, "Calculate")
        botao(left, exemplo_derivada_implicita, "Example")

        resultado_text_implicita = ctk.CTkTextbox(right, font=font)
        resultado_text_implicita.pack(fill="both", expand=True)

    # ====================== ABA LHOPITAL =========================
    def aba_lhopital(self, frame):
        global entrada_num, entrada_den, entrada_ponto, direcao_lhopital, resultado_text_lhopital

        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("When to use L'Hôpital?"), command=abrir_explicacao_lhopital).pack(pady=5, anchor="w")

        entrada_num = labeled_input(left, "Numerator Function:")
        aplicar_validacao_em_tempo_real(entrada_num)
        entrada_den = labeled_input(left, "Denominator Function:")
        aplicar_validacao_em_tempo_real(entrada_den)
        entrada_ponto = labeled_input(left, "Approaching:")
        aplicar_validacao_em_tempo_real(entrada_ponto)

        direcao_lhopital = ctk.StringVar(value="+")
        ctk.CTkOptionMenu(left, variable=direcao_lhopital, values=["+", "-"]).pack(pady=5, anchor="w")

        botao(left, calculo_lhopital, "Apply L'Hôpital")
        botao(left, exemplo_lhopital, "Example")

        resultado_text_lhopital = ctk.CTkTextbox(right, font=font)
        resultado_text_lhopital.pack(fill="both", expand=True)


    # ====================== ABA INTEGRAIS =========================
    def aba_integrais(self, frame):
        global entrada_integrais, entrada_limite_inf, entrada_limite_sup, resultado_text_integral
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text=_("What is an Integral?"), command=abrir_explicacao_integral).pack(pady=5, anchor="w")

        entrada_integrais = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entrada_integrais)
        entrada_limite_inf = labeled_input(left, "Lower limit:")
        aplicar_validacao_em_tempo_real(entrada_limite_inf)
        entrada_limite_sup = labeled_input(left, "Upper limit:")
        aplicar_validacao_em_tempo_real(entrada_limite_sup)

        botao(left, calculo_integral, "Calculate")
        botao(left, exemplo_integral, "Example")

        resultado_text_integral = ctk.CTkTextbox(right, font=font)
        resultado_text_integral.pack(fill="both", expand=True)

        
        try:
            nome_imagem = "integral2.png" if idioma_atual == "PT" else "integral.png"
            img = ctk.CTkImage(Image.open(nome_imagem), size=(300, 180))
            ctk.CTkLabel(right, image=img, text="").pack(pady=10)
        except Exception:
            pass

    # ====================== ABA MANUAL =========================
    def aba_manual(self, frame):
        ctk.CTkButton(
            frame,
            text=_("Open DDX Manual"),
            command=lambda: webbrowser.open('https://docs.google.com/document/d/1hvcUL36juGBm_8lsdOpPrMLWzmYnGvakKHaMj1BbxlY/edit?usp=sharing'),
            width=300
        ).pack(pady=20)


# ====================== EXECUÇÃO =========================
if __name__ == "__main__":
    initial_page = InitialPage()
    initial_page.mainloop()