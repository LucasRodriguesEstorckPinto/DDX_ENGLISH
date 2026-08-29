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


dados_x = None
dados_y = None
interpolar_var = None
botao_plot_dados = None
check_interpolar = None


matplotlib.use("TkAgg")

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("green")

# Fonte padrão para os widgets
font = ("Segoe UI", 14)

# Configurações globais
sp.init_printing()
x = sp.symbols('x')
n = sp.symbols('n', integer=True)


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
                f"The {i}-th order derivative of the function is: {derivada_atual}\n"
            )

        
        point_str = entradaponto.get()
        if point_str:
            point = float(sp.sympify(point_str))
            valor_derivada = derivada_atual.subs(x, point)

            resultado_text_deriv.insert(
                ctk.END,
                f"\nThe value of the {ordem}-th derivative at point x={point} is: {valor_derivada}\n"
            )

            
            if ordem == 1:
                coef_angular = valor_derivada
                reta = func.subs(x, point) + coef_angular * (x - point)
                resultado_text_deriv.insert(
                    ctk.END,
                    f"The equation of the tangent line is: {reta}\n\n"
                )

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating the derivative: {e}")

def calculo_derivada_implicita():
    global entrada_implicita, entrada_var_dependente, entrada_var_independente, resultado_text_implicita
    try:
        
        eq_str = entrada_implicita.get()
        y_str = entrada_var_dependente.get().strip()
        x_str = entrada_var_independente.get().strip()

        # Validação básica
        if not eq_str or not y_str or not x_str:
            messagebox.showerror("Error", "Please fill in all fields: equation, dependent variable (y) and independent variable (x).")
            return

        
        x_sym = sp.symbols(x_str)
        y_sym = sp.symbols(y_str)

        
        if '=' in eq_str:
            lhs_str, rhs_str = eq_str.split('=', 1)
            lhs = sp.sympify(lhs_str.strip())
            rhs = sp.sympify(rhs_str.strip())
            eq = sp.Eq(lhs, rhs)
        else:
            
            eq = sp.sympify(eq_str)

       
        derivada_implicita = sp.idiff(eq, y_sym, x_sym)

        
        resultado_text_implicita.delete("1.0", ctk.END)

        
        resultado_text_implicita.insert(ctk.END, f"The derivative of {y_str} with respect to {x_str} (d{y_str}/d{x_str}) is:\n\n")
        resultado_text_implicita.insert(ctk.END, f"{derivada_implicita}")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating the implicit derivative.\nCheck the equation and the variables.\n\nDetails: {e}")


def calculo_limite():
    global resultado_text_limite, entradalimit, entradavar, entradatend, direcao_var
    try:
        func_str = entradalimit.get()
        func = sp.sympify(func_str)
        variavel = sp.symbols(entradavar.get())
        valor_tendencia = float(sp.sympify(entradatend.get()))
        direcao = direcao_var.get()  

        if direcao == "Both":
            limite_esquerda = sp.limit(func, variavel, valor_tendencia, dir='-')
            limite_direita = sp.limit(func, variavel, valor_tendencia, dir='+')

            resultado_text_limite.delete("1.0", ctk.END)
            if limite_esquerda == limite_direita:
                resultado_text_limite.insert(ctk.END, f"The limit of the function is: {limite_esquerda}")
            else:
                resultado_text_limite.insert(ctk.END, f"The limit of the function does not exist.")
        else:
            if direcao == "Left":
                limite = sp.limit(func, variavel, valor_tendencia, dir='-')
            elif direcao == "Right":
                limite = sp.limit(func, variavel, valor_tendencia, dir='+')

            resultado_text_limite.delete("1.0", ctk.END)
            resultado_text_limite.insert(ctk.END, f"The limit of the function from the {direcao.lower()} is: {limite}")

    except Exception as e:
        messagebox.showerror("Error", "An error occurred while calculating the limit. Check your input.")

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
        resultado_text_raiz.insert(ctk.END, f"The {indice}-th root of {numero} is: {raiz_value:.4}")
    except ValueError:
        messagebox.showerror("Error", "Please provide a valid index and/or number to calculate the root.")



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
                linha = re.sub(r'[,\t]+', ' ', linha.strip())  # trata vírgulas e tabs
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
        resultado_text_grafico.insert(ctk.END, "\nFile loaded successfully. Points ready for plotting.\n")

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        messagebox.showerror("Error importing file", f"Error processing file: {str(e)}")

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
                    resultado_text_grafico.insert(ctk.END, "\nInterpolated curve plotted successfully.\n")
                except Exception as e:
                    import traceback
                    print(traceback.format_exc())
                    messagebox.showwarning("Interpolation failed", f"Error interpolating: {str(e)}")
            else:
                resultado_text_grafico.insert(ctk.END, "\n⚠️ At least 4 points are required for cubic interpolation.\n")

        
        ax.scatter(x_ord, y_ord, color='red', s=60, zorder=5, label="Points")

    
        ax.axhline(0, color='black', lw=1.2, linestyle='dashed')
        ax.axvline(0, color='black', lw=1.2, linestyle='dashed')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title("Imported data plot")
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

        resultado_text_grafico.insert(ctk.END, "\nPlot generated.\n")

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        messagebox.showerror("Error", f"Error plotting points: {str(e)}")



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
                raise ValueError("'Piecewise' mode selected. Insert functions and intervals.")
            if len(func_list_str) != len(interval_list_str):
                raise ValueError("The number of functions and intervals (separated by ';') must be the same.")

            for f_str, i_str in zip(func_list_str, interval_list_str):
                func_sym = sp.sympify(f_str)

              
                parts = [p.strip() for p in i_str.split(',')]
                if len(parts) != 2:
                    raise ValueError(f"Badly formatted interval: '{i_str}'. Use 'a, b'.")

                lower = float(sp.N(sp.sympify(parts[0])))
                upper = float(sp.N(sp.sympify(parts[1])))
                if lower > upper:
                    lower, upper = upper, lower

                lista_de_pedacos.append( (func_sym, lower, upper) )
                full_lower = min(full_lower, lower)
                full_upper = max(full_upper, upper)

        else:
            
            func_list_str = [f.strip() for f in func_input_str.split(',') if f.strip()]
            if not func_list_str:
                raise ValueError("No function inserted.")

           
            parts = [p.strip() for p in intervalo_input_str.split(',')]
            if len(parts) != 2:
                raise ValueError(f"Badly formatted interval: '{intervalo_input_str}'. Use 'a, b'.")

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

        if y_lower is not None:
             result_text = f'Y interval used: [{y_lower:.2f}, {y_upper:.2f}]\n' + result_text


        

        for i, (func_sym, lower, upper) in enumerate(lista_de_pedacos):

           
            result_text += f"\n--- Analysis of ${sp.latex(func_sym)}$ in [{lower}, {upper}] ---\n"

          
            func_numeric = sp.lambdify(x, func_sym, 'numpy')
            x_vals = ajustar_amostragem(lower, upper) 

           
            y_vals = func_numeric(x_vals)
            ax.plot(x_vals, y_vals, label=f'${sp.latex(func_sym)}$ in [{lower},{upper}]', linewidth=2.5, color=f'C{i}')

           
            try:
                if func_sym.has(sp.tan):
                    n_vals = range(int(lower/sp.pi)-1, int(upper/sp.pi)+2)
                    vertical_asymptotes = [sp.pi/2 + n*sp.pi for n in n_vals]
                else:
                    vertical_asymptotes = sp.singularities(func_sym, x)
                vertical_asymptotes = [asy for asy in vertical_asymptotes if asy.is_real]
                for asy in vertical_asymptotes:
                    asy_val = float(asy.evalf())
                    if lower < asy_val < upper: 
                        ax.axvline(asy_val, color='magenta', linestyle='--', linewidth=2)
                        result_text += f'Vertical asymptote at x = {asy_val:.2f}\n'
            except Exception as e:
                print(f"Error calculating vertical asymptotes: {e}")

            
            try:
                lim_neg = sp.limit(func_sym, x, -sp.oo)
                lim_pos = sp.limit(func_sym, x, sp.oo)
                for lim, side in [(lim_neg, '-∞'), (lim_pos, '+∞')]:
                    if getattr(lim, 'is_real', False) and not lim.has(sp.oo, sp.zoo):
                        lim_val = float(lim.evalf())
                        ax.axhline(lim_val, color='cyan', linestyle='--', linewidth=2)
                        result_text += f'Horizontal asymptote at y = {lim_val:.2f} (limit at {side})\n'
            except Exception as e:
                print(f"Error calculating horizontal asymptotes: {e}")

            
            try:
                coef, intercept = encontrar_assintota_obliqua(func_sym, x)
                if coef is not None and intercept is not None:
                    ax.axline((0, float(intercept)), slope=float(coef), color='orange', linestyle='--')
                    result_text += f'Oblique asymptote: y = {float(coef):.2f}x + {float(intercept):.2f}\n'
            except Exception as e:
                print(f"Error calculating oblique asymptote: {e}")

            
            try:
                fprime, fsecond = calcular_derivadas(func_sym, x)
                
                cp = numerical_roots(fprime, x, lower, upper)
                ip = numerical_roots(fsecond, x, lower, upper)

                if show_points_var.get():
                    colors = ['#e41a1c', '#4daf4a', '#ff7f00', '#984ea3', '#377eb8']
                    markers = ['^', 'v', 'D', 'o', 's']
                    for p, color, marker in zip(cp, colors[:len(cp)], markers[:len(cp)]):
                        try:
                            y_p = float(func_sym.subs(x, p).evalf())
                            fsecond_val = float(fsecond.subs(x, p).evalf())
                            point_type = "Maximum" if fsecond_val < 0 else "Minimum" if fsecond_val > 0 else "Saddle"
                        except Exception:
                            y_p = float(func_sym.subs(x, p).evalf())
                            point_type, color, marker = "Critical", '#984ea3', 'o'

                        ax.scatter(p, y_p, color=color, marker=marker, s=100, edgecolors='black', zorder=6)
                        ax.annotate(
                            f'{point_type}\n({p:.2f}, {y_p:.2f})',
                            xy=(p, y_p), xytext=(p + 0.4, y_p + (0.4 if point_type == "Maximum" else -0.4)),
                            textcoords='data', fontsize=10, fontweight='bold', color='white',
                            bbox=dict(boxstyle='round,pad=0.3', fc=color, ec='none'),
                            arrowprops=dict(arrowstyle='-|>', color=color, lw=1.5), zorder=7
                        )
                        result_text += f'Local {point_type.lower()} at ({p:.2f}, {y_p:.2f})\n'

                    for p in ip:
                        y_p = float(func_sym.subs(x, p).evalf())
                        ax.scatter(p, y_p, color='#377eb8', marker='s', s=100, edgecolors='black', zorder=6)
                        ax.annotate(
                            f'Inflection\n({p:.2f}, {y_p:.2f})',
                            xy=(p, y_p), xytext=(p + 0.4, y_p + 0.4), textcoords='data',
                            fontsize=10, fontweight='bold', color='white',
                            bbox=dict(boxstyle='round,pad=0.3', fc='#377eb8', ec='none'),
                            arrowprops=dict(arrowstyle='-|>', color='#377eb8', lw=1.5), zorder=7
                        )
                        result_text += f'Local inflection at ({p:.2f}, {y_p:.2f})\n'
                else:
                    if i == 0: 
                        result_text += "Points not explicitly shown (checkbox disabled).\n"

               
                growth_points = sorted(list(set([lower] + cp + [upper])))
                for j in range(len(growth_points) - 1):
                    mid = (growth_points[j] + growth_points[j+1]) / 2
                    try:
                        derivative_mid = float(fprime.subs(x, mid).evalf())
                        if derivative_mid > 0:
                            result_text += f'Growth in [{growth_points[j]:.2f}, {growth_points[j+1]:.2f}]\n'
                        elif derivative_mid < 0:
                            result_text += f'Decrease in [{growth_points[j]:.2f}, {growth_points[j+1]:.2f}]\n'
                        else:
                            result_text += f'Constant in [{growth_points[j]:.2f}, {growth_points[j+1]:.2f}]\n'
                    except Exception:
                        continue

                
                if show_points_var.get():
                    try:
                        conc_points = [float(p) for p in ip]
                        sing2 = [s for s in sp.singularities(fsecond, x) if s.is_real]
                        sing2 = [float(s.evalf()) for s in sing2]

                        internal_points = sorted(set([p for p in conc_points + sing2 if lower < p < upper]))
                       
                        breakpoints = sorted(list(set([lower] + internal_points + [upper])))

                        conc_up_labeled = False
                        conc_down_labeled = False

                        for j in range(len(breakpoints) - 1):
                            a = breakpoints[j]
                            b = breakpoints[j + 1]
                            if a == b: continue
                            mid = (a + b) / 2

                            try:
                                val = float(fsecond.subs(x, mid).evalf())
                            except Exception:
                                try:
                                    f2_numeric = sp.lambdify(x, fsecond, 'numpy')
                                    val = float(f2_numeric(mid))
                                except Exception:
                                    continue

                            if val > 0:
                                ax.axvspan(a, b, alpha=0.12, facecolor='green', zorder=1,
                                           label='Positive concavity' if not conc_up_labeled else None)
                                conc_up_labeled = True
                                result_text += f'Positive concavity in [{a:.2f}, {b:.2f}]\n'
                            elif val < 0:
                                ax.axvspan(a, b, alpha=0.12, facecolor='red', zorder=1,
                                           label='Negative concavity' if not conc_down_labeled else None)
                                conc_down_labeled = True
                                result_text += f'Negative concavity in [{a:.2f}, {b:.2f}]\n'
                    except Exception as e:
                        print(f"Error calculating concavity intervals: {e}")

            except Exception as e:
                result_text += f"Error analyzing derivatives for this piece: {e}\n"

            
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

       
        resultado_text_grafico.delete("1.0", ctk.END)
        resultado_text_grafico.insert(ctk.END, result_text + "\nPlot generated successfully!")
        resultado_text_grafico.insert(ctk.END, "\nNote: All points shown are LOCAL to their respective interval.")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while plotting the graph: {str(e)}")


def validar_entrada(func_str):
    pattern = r'^[a-zA-Z0-9\s+\-*/().^sincoslogexp]+$'
    if not re.match(pattern, func_str):
        raise ValueError("Invalid input: use only valid mathematical characters.")
    return func_str.replace("^", "**").replace("sen", "sin").replace("arctg", "atan").replace("arcsen", "asin").replace("arccos", "acos")

def formatar_intervalo(intervalo):
    if not isinstance(intervalo, (sp.Interval, sp.Union, str)):
        return str(intervalo)

    if isinstance(intervalo, str):
        return intervalo

    if isinstance(intervalo, sp.Union):
        intervalos_formatados = [formatar_intervalo(i) for i in intervalo.args]
        return " ∪ ".join(intervalos_formatados)

    esquerda = intervalo.left
    direita = intervalo.right

    if esquerda == -sp.oo:
        inicio = "(-∞"
    else:
        valor_esq = float(esquerda.evalf())
        if abs(valor_esq - round(valor_esq)) < 1e-10:
            valor_esq = int(round(valor_esq))
        elif abs(valor_esq) < 1000:
            valor_esq_str = str(round(valor_esq, 4))
            valor_esq = valor_esq_str.rstrip('0').rstrip('.') if '.' in valor_esq_str else valor_esq_str

        inicio = f"[{valor_esq}" if not intervalo.left_open else f"({valor_esq}"

    if direita == sp.oo:
        fim = "+∞)"
    else:
        valor_dir = float(direita.evalf())
        if abs(valor_dir - round(valor_dir)) < 1e-10:
            valor_dir = int(round(valor_dir))
        elif abs(valor_dir) < 1000:
            valor_dir_str = str(round(valor_dir, 4))
            valor_dir = valor_dir_str.rstrip('0').rstrip('.') if '.' in valor_dir_str else valor_dir_str

        fim = f"{valor_dir}]" if not intervalo.right_open else f"{valor_dir})"

    return f"{inicio}, {fim}"


def formatar_conjunto(conjunto):
    
    if isinstance(conjunto, str):
        return conjunto
    try:
        return str(sp.simplify(conjunto))
    except Exception:
        return repr(conjunto)


def explicar_dominio(dominio, func_str=""):

    if isinstance(dominio, str):
        return dominio  

  
    if dominio == S.Reals:
        return "All real numbers (no restrictions on the function)."

    explicacao = []

    if isinstance(dominio, Complement):
        conj_principal = dominio.args[0]
        conj_excluido = dominio.args[1]
        if isinstance(conj_excluido, FiniteSet):
            pontos = ", ".join([f"x ≠ {p}" for p in conj_excluido])
            explicacao.append(f"The function is not defined at {pontos} due to the presence of division by zero.")
        dominio = conj_principal

    if isinstance(dominio, Union):
        for intervalo in dominio.args:
            explicacao.append(f"Allowed interval: {intervalo}")
    elif isinstance(dominio, Interval):
        if dominio.left == 0 and dominio.right == oo:
            explicacao.append("The domain is (0, +∞) because the argument of the logarithm or denominator must be positive.")
        elif dominio.left == 1 and dominio.right == oo:
            explicacao.append("The domain is (1, +∞) because the argument of ln(x - 1) must be greater than zero.")
        elif dominio.left == 3 and dominio.right == oo:
            explicacao.append("The domain is [3, +∞) because the square root requires the expression inside it to be non-negative.")
        else:
            explicacao.append(f"The domain is {dominio} due to symbolic restrictions of the function.")

    elif isinstance(dominio, FiniteSet):
        explicacao.append(f"Domain restricted to specific values: {', '.join(str(v) for v in dominio)}")

    return " ".join(explicacao)


def explicar_imagem(imagem, func_str):
    if isinstance(imagem, str):
        return imagem

    if "sin" in func_str or "cos" in func_str:
        if imagem == sp.Interval(-1, 1):
            return "The range is [-1, 1]. Sine and cosine functions vary only between -1 and 1."

    if "tan" in func_str:
        if "ℝ" in str(imagem):
            return "The range is ℝ (all real numbers). The tangent function can take any real value."

    try:
        if "x**" in func_str:
            if "x**2" in func_str and imagem == sp.Interval(0, sp.oo):
                return "The range is [0, +∞). Quadratic functions of the form ax² + bx + c (with a > 0) are bounded below."
            if imagem == sp.S.Reals:
                return "The range is ℝ (all real numbers). This polynomial function takes all real values."
    except Exception:
        pass

    return f"The range of the function is {formatar_conjunto(imagem)}."


def calcular_dominio(func, x):
   
    try:
        func = sp.sympify(func)
        dominio = S.Reals
        restricoes = []

        try:
            sins = sp.singularities(func, x)
            
            if isinstance(sins, (sp.FiniteSet, set, list, tuple)):
                for s in sins:
                    if getattr(s, 'is_real', False):
                        dominio = dominio - FiniteSet(sp.nsimplify(s))
        except Exception:
           
            pass

       
        try:
            _, denom = func.as_numer_denom()
            if denom != 1:
                zeros = sp.solveset(denom, x, domain=sp.S.Reals)
                if isinstance(zeros, sp.FiniteSet):
                    dominio = dominio - zeros
                else:
                    
                    pass
        except Exception:
            pass

       
        for p in func.atoms(sp.Pow):
            exp = p.exp
            base = p.base
            
            if getattr(exp, 'is_Rational', False) and (exp.q % 2 == 0):
               
                try:
                    sol = solve_univariate_inequality(base >= 0, x)
                    restricoes.append(sol)
                except Exception:
                    pass

        
        for l in func.atoms(sp.log):
            arg = l.args[0]
            try:
                sol = solve_univariate_inequality(arg > 0, x)
                restricoes.append(sol)
            except Exception:
                pass

        
        for r in restricoes:
            try:
                dominio = dominio.intersect(r)
            except Exception:
                
                pass

        
        try:
            dominio = sp.simplify(dominio)
        except Exception:
            pass

        return dominio

    except Exception as e:
        return f"Error calculating the domain: {e}"

def calcular_imagem(func, x, dominio):
   
    try:
        func = sp.simplify(sp.sympify(func))

        
        if isinstance(dominio, str):
            return f"Could not calculate range because the domain is invalid: {dominio}"

        
        if func.has(sp.Abs):
            
            for atom in func.atoms(sp.Abs):
                inner = sp.simplify(atom.args[0])
                
                try:
                    poly_inner = sp.Poly(inner, x)
                    if poly_inner.degree() == 1:
                        a = float(poly_inner.coeffs()[0])  # coef de x
                        b = float(poly_inner.coeffs()[1]) if len(poly_inner.coeffs()) > 1 else 0.0
                       
                        if abs(a) > 0:
                            root = -b / a
                            
                            try:
                                if hasattr(dominio, 'contains') and dominio.contains(sp.nsimplify(root)):
                                    
                                    return sp.Interval(0, sp.oo)
                                else:
                                    
                                    pass
                            except Exception:
                               
                                return sp.Interval(0, sp.oo)
                except Exception:
                    pass
            
        try:
            rng = function_range(func, x, domain=dominio)
            if not isinstance(rng, sp.ConditionSet):
                return sp.simplify(rng)
        except Exception:
            rng = None

        # helpers para infinito
        def _is_pos_inf(v):
            return v == sp.oo or v == sp.zoo
        def _is_neg_inf(v):
            return v == -sp.oo or v == sp.zoo

        
        unbounded_right = False
        unbounded_left = False
        try:
            lim_plus = sp.limit(func, x, sp.oo)
            lim_minus = sp.limit(func, x, -sp.oo)
            if _is_pos_inf(lim_plus) or _is_pos_inf(lim_minus):
                unbounded_right = True
            if _is_neg_inf(lim_plus) or _is_neg_inf(lim_minus):
                unbounded_left = True
        except Exception:
            pass

        try:
            sins = sp.singularities(func, x)
            for s in (list(sins) if isinstance(sins, (list, tuple, set)) else ([sins] if sins else [])):
                if getattr(s, 'is_real', False):
                    try:
                        in_dom = True
                        if hasattr(dominio, 'contains'):
                            in_dom = dominio.contains(s)
                    except Exception:
                        in_dom = True
                    if in_dom:
                        try:
                            l_plus = sp.limit(func, x, s, dir='+')
                            l_minus = sp.limit(func, x, s, dir='-')
                            if _is_pos_inf(l_plus) or _is_pos_inf(l_minus):
                                unbounded_right = True
                            if _is_neg_inf(l_plus) or _is_neg_inf(l_minus):
                                unbounded_left = True
                        except Exception:
                            unbounded_right = True
                            unbounded_left = True
        except Exception:
            pass

       
        sample_y = []

        
        intervals = []
        points_to_force = []

        if isinstance(dominio, sp.Interval):
            intervals = [dominio]
        elif isinstance(dominio, sp.Union):
            intervals = [i for i in dominio.args]
        elif isinstance(dominio, sp.FiniteSet):
           
            pts = []
            for p in dominio:
                try:
                    pts.append(float(p.evalf()))
                except Exception:
                    pass
            for p in pts:
                try:
                    yv = func.subs(x, p).evalf()
                    yvf = float(yv)
                    if np.isfinite(yvf):
                        sample_y.append(yvf)
                except Exception:
                    pass
            intervals = []
        else:
            intervals = [sp.Interval(-100.0, 100.0)]

        
        try:
            if hasattr(dominio, 'contains') and dominio.contains(0):
                points_to_force.append(0.0)
        except Exception:
           
            points_to_force.append(0.0)

        
        try:
            deriv = sp.diff(func, x)
            cps = sp.solve(sp.Eq(deriv, 0), x)
            for cp in cps:
                try:
                    if float(sp.re(cp)) == float(cp):  # real-ish
                        points_to_force.append(float(sp.N(cp)))
                except Exception:
                    pass
        except Exception:
            pass

        
        BIG = 1000.0
        for interval in intervals:
            try:
                left = interval.left
                right = interval.right
            except Exception:
                left, right = -100.0, 100.0

            a = float(left) if left != -sp.oo else -min(BIG, 1000.0)
            b = float(right) if right != sp.oo else min(BIG, 1000.0)

            if a == b:
                xs = np.array([a])
            else:
                span = b - a
                n = int(min(max(300, int(abs(span) * 80)), 3000))
                xs = np.linspace(a, b, n)

            
            forced_in_interval = [p for p in points_to_force if a - 1e-12 <= p <= b + 1e-12]
            if forced_in_interval:
                xs = np.unique(np.concatenate((xs, np.array(forced_in_interval))))

            for xv in xs:
                try:
                    
                    yv = func.subs(x, float(xv))
                    yv_num = float(sp.N(yv))
                    if np.isfinite(yv_num):
                        sample_y.append(yv_num)
                    else:
                        if yv == sp.oo or yv == sp.zoo:
                            unbounded_right = True
                        if yv == -sp.oo:
                            unbounded_left = True
                except Exception:
                    
                    continue

        
        for p in points_to_force:
            try:
                yv = func.subs(x, float(p))
                yv_num = float(sp.N(yv))
                if np.isfinite(yv_num):
                    sample_y.append(yv_num)
            except Exception:
                pass

        if not sample_y:
            return "Undefined range (no valid values found by sampling)"

        min_val = min(sample_y)
        max_val = max(sample_y)

        
        try:
            zero_in_domain = hasattr(dominio, 'contains') and dominio.contains(0)
        except Exception:
            zero_in_domain = False

        TOL = 1e-9
        if zero_in_domain and abs(min_val) <= TOL:
            min_val = 0.0

       
        if unbounded_left and unbounded_right:
            return sp.S.Reals
        if unbounded_right:
            left_val = sp.nsimplify(min_val) if abs(min_val) > 0 else sp.Integer(0)
            
            return sp.Interval(left_val, sp.oo)
        if unbounded_left:
            right_val = sp.nsimplify(max_val) if abs(max_val) > 0 else sp.Integer(0)
            return sp.Interval(-sp.oo, right_val)

        
        if abs(max_val - min_val) < 1e-12:
            return FiniteSet(sp.nsimplify(min_val))
        return sp.Interval(min_val, max_val)

    except Exception as e:
        return f"Error calculating the range: {e}"


def calculo_dominio_imagem():
   
    global resultado_text_dom, entradadom, grafico_label, x

    try:
        func_str = entradadom.get()
        func_str = validar_entrada(func_str)  

        func = sp.sympify(func_str)

       
        dominio = calcular_dominio(func, x)

        if isinstance(dominio, str) and "Error" in dominio:
            imagem = f"Could not calculate range because the domain is invalid: {dominio}"
        else:
            imagem = calcular_imagem(func, x, dominio)

        
        dominio_explicado = explicar_dominio(dominio, func_str)
        imagem_explicada = explicar_imagem(imagem, func_str)

        resultado = f"""Results:
========================
Function: {func_str}

Domain: {formatar_conjunto(dominio)}
{dominio_explicado}

Range: {formatar_conjunto(imagem)}
{imagem_explicada}
========================"""

        resultado_text_dom.delete("1.0", ctk.END)
        resultado_text_dom.insert(ctk.END, resultado)

        
        if 'grafico_label' in globals() and grafico_label is not None:
            try:
                grafico_label.destroy()
            except Exception:
                pass

        
        try:
            x_vals = np.linspace(-10, 10, 1000)
            x_plot = []
            y_plot = []
            for val in x_vals:
                try:
                    yv = func.subs(x, float(val)).evalf()
                    yvf = float(yv)
                    if np.isfinite(yvf):
                        x_plot.append(val)
                        y_plot.append(yvf)
                except Exception:
                    continue

            if x_plot and y_plot:
                plt.figure(figsize=(10, 6))
                plt.plot(x_plot, y_plot, label=func_str)
                plt.title(f"Graph of {func_str}")
                plt.xlabel("x")
                plt.ylabel("f(x)")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.savefig("grafico.png")
                plt.close()

                img = ctk.CTkImage(Image.open("grafico.png"), size=(400, 300))
                grafico_label = ctk.CTkLabel(master=resultado_text_dom.master, image=img, text="")
                grafico_label.pack(pady=10, after=resultado_text_dom)
            else:
               
                print("No valid points for plotting (restricted domain or many evaluation errors).")

        except Exception as e:
            print(f"Error generating graph: {e}")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating domain and range: {str(e)}")

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
            resultado_text_integral.insert(ctk.END, f"The definite integral of the function from {limite_inf} to {limite_sup} is: {integral_def}\n")
        else:
            integral = sp.integrate(func, x)
            resultado_text_integral.delete("1.0", ctk.END)
            resultado_text_integral.insert(ctk.END, f"The indefinite integral of the function is: {integral} + C\n")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred while calculating the integral:")

def calculo_derivada_implicita():
        global entrada_implicita, entrada_var_dependente, entrada_var_independente, resultado_text_implicita
        try:
            
            eq_str = entrada_implicita.get()
            y_str = entrada_var_dependente.get().strip()
            x_str = entrada_var_independente.get().strip()

            
            if not eq_str or not y_str or not x_str:
                messagebox.showerror("Error", "Please fill in all fields: equation, dependent variable (y) and independent variable (x).")
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

           
            resultado_text_implicita.insert(ctk.END, f"The derivative of {y_str} with respect to {x_str} (d{y_str}/d{x_str}) is:\n\n")
            resultado_text_implicita.insert(ctk.END, f"{derivada_implicita}")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while calculating the implicit derivative.\nCheck the equation and the variables.\n\nDetails: {e}")


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
        messagebox.showerror("Error", f"An error occurred while plotting the graph:")

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


def exemplo_lhopital(self):
    entrada_num.delete(0, ctk.END)
    entrada_den.delete(0, ctk.END)
    entrada_ponto.delete(0, ctk.END)
    entrada_num.insert(0, "sin(x)")
    entrada_den.insert(0, "x")
    entrada_ponto.insert(0, "0")

def exemplo_raiz():
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
    example_text = (
        "Example of Integral:\n"
        "Function: f(x) = x^2\n"
        "Indefinite Integral: ∫x^2 dx = (1/3)x^3 + C, where C is the constant of integration.\n"
        "Definite Integral from 0 to 2: ∫(from 0 to 2) x^2 dx = [(1/3)x^3] from 0 to 2 = (8/3) - 0 = 8/3.\n"
        "This represents the area under the curve of f(x) between x=0 and x=2."
    )
    resultado_text_integral.delete("1.0", ctk.END)
    resultado_text_integral.insert(ctk.END, example_text)

def abrir_explicacao_integral():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title("Explanation of Integrals")
    janela_explicacao.geometry("500x300")

    texto_explicacao = """The integral of a function represents the area under the curve of that function in a given interval.
It is used to calculate areas, volumes, and solve physical problems such as work and displacement.

Source: Stewart, James. Calculus. 8th edition."""

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text="Close", command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_derivada():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title("Explanation of Derivatives")
    janela_explicacao.geometry("500x300")

    texto_explicacao = """The derivative of a function represents the rate of change of that function at a given point.
It is used to calculate velocities, accelerations, and solve physical problems such as optimization and population growth.

Source: Stewart, James. Calculus. 8th edition."""

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text="Close", command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_limites():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title("Explanation of Limits")
    janela_explicacao.geometry("500x300")

    texto_explicacao = """The limit of a function describes the behavior of that function as the independent variable approaches a certain value.
It is used to define derivatives, integrals, and solve problems involving continuity and asymptotic behavior.

Source: Stewart, James. Calculus. 8th edition."""

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text="Close", command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_derivadas_parciais():
    pass

def abrir_explicacao_dominios():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title("Explanation of Domains")
    janela_explicacao.geometry("500x300")

    texto_explicacao = """The domain of a function is the set of all input values for which the function is defined.
The range of a function is the set of all output values that the function can take.
They are used to understand the behavior and restrictions of functions in various mathematical and applied contexts.

Source: Stewart, James. Calculus. 8th edition."""

    label_texto = ctk.CTkLabel(janela_explicacao, text=texto_explicacao, wraplength=450, justify="left")
    label_texto.pack(padx=20, pady=20)

    botao_fechar = ctk.CTkButton(janela_explicacao, text="Close", command=janela_explicacao.destroy)
    botao_fechar.pack(pady=10)

def abrir_explicacao_lhopital():
    janela_explicacao = ctk.CTkToplevel()
    janela_explicacao.title("Explanation of L'Hôpital's Rule")
    janela_explicacao.geometry("600x320")

    texto = """L'Hôpital's Rule is used to solve limits that present indeterminate forms,
such as 0/0 or ∞/∞.

Let f(x) and g(x) be differentiable functions on an open interval containing 'a', and if:
- lim(x→a) f(x) = 0 and lim(x→a) g(x) = 0  or
- lim(x→a) f(x) = ∞ and lim(x→a) g(x) = ∞

Then:
    lim(x→a) f(x)/g(x) = lim(x→a) f'(x)/g'(x)
provided this derivative limit exists.

The rule can be applied repeatedly until the indeterminacy disappears."""

    label = ctk.CTkLabel(janela_explicacao, text=texto, justify="left", wraplength=580, font=("Segoe UI", 14))
    label.pack(padx=20, pady=20)

class ModernEntry(ctk.CTkEntry):
    def get(self):
        text = super().get()
       
        text = re.sub(r'(?<=\d)(?=pi\b)', '*', text, flags=re.IGNORECASE)
        text = re.sub(r'(?<=\d)(?=e\b)', '*', text, flags=re.IGNORECASE)
        
        text = re.sub(r'\bpi\b', 'pi', text, flags=re.IGNORECASE)
        text = re.sub(r'\be\b', 'E', text, flags=re.IGNORECASE)
        return text

# Cria um rótulo e uma entrada logo abaixo
def labeled_input(parent, label_text):
    frame = ctk.CTkFrame(parent)
    frame.pack(anchor="w", pady=5, padx=5, fill="x")

    label = ctk.CTkLabel(frame, text=label_text, font=font)
    label.pack(anchor="w", padx=5)

    entry = ctk.CTkEntry(frame, width=400, height=30, corner_radius=5, font=font)
    entry.pack(padx=5, pady=5, anchor="w")

    return entry


def botao(parent, func, texto):
    btn = ctk.CTkButton(parent, text=texto, command=func, width=200, corner_radius=5, font=font)
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
        self.geometry("500x300")
        self.resizable(False, False)
        self.configure(padx=20, pady=20)

        ctk.CTkLabel(self, text="Welcome to the DDX Calculator", font=("Segoe UI", 20, "bold")).pack(pady=20)

        open_calculator_btn = ctk.CTkButton(self, text="Open DDX Calculator", command=self.open_calculator, width=250)
        open_calculator_btn.pack(pady=10)

        manual_btn = ctk.CTkButton(
            self,
            text="Open DDX Manual",
            command=lambda: webbrowser.open('https://drive.google.com/file/d/1XhUZMxmc4bNwYTh5FOoXdG_eztrMHVSs/view?usp=sharing'),
            width=250
        )
        manual_btn.pack(pady=10)

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

        abas = ["Domain and Range", "Limits", "Derivatives", "Root", "Graphs", "Implicit Derivatives", "L'Hôpital", "Integrals", "Manual"]
        frames = {aba: tabview.add(aba) for aba in abas}

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

        ctk.CTkButton(left, text="What are Domains and Ranges?", command=abrir_explicacao_dominios).pack(pady=5, anchor="w")

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

        ctk.CTkButton(left, text="What is a Derivative?", command=abrir_explicacao_derivada).pack(pady=5, anchor="w")

        entradaderiv = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entradaderiv)
        entradaponto = labeled_input(left, "Point:")
        entradaordem = labeled_input(left, "Derivative order (e.g., 1, 2, 3...):")

        botao(left, calculo_derivada, "Calculate")
        botao(left, exemplo_derivada, "Example")
        botao(left, plot_func_tangente, "Plot Tangent (Order 1)")

        resultado_text_deriv = ctk.CTkTextbox(right, font=font)
        resultado_text_deriv.pack(fill="both", expand=True)

        img = ctk.CTkImage(Image.open("deriva.png"), size=(250, 120))
        ctk.CTkLabel(right, image=img, text="").pack(pady=10)


    # ====================== ABA DERIVADAS PARCIAIS =========================
    def aba_derivadas_parciais(self, frame):
        global entradafuncparcial, entradavarparcial, resultado_text_parcial
        left, right = self.estrutura_aba(frame)

        ctk.CTkButton(left, text="What are Partial Derivatives?", command=abrir_explicacao_derivadas_parciais).pack(pady=5, anchor="w")

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

        ctk.CTkButton(left, text="What are Limits?", command=abrir_explicacao_limites).pack(pady=5, anchor="w")

        entradalimit = labeled_input(left, "Function:")
        aplicar_validacao_em_tempo_real(entradalimit)
        entradavar = labeled_input(left, "Variable:")
        entradatend = labeled_input(left, "Approaching:")
        aplicar_validacao_em_tempo_real(entradatend)

        direcao_var = ctk.StringVar(value="Both")
        ctk.CTkOptionMenu(left, variable=direcao_var, values=["Left", "Right", "Both"]).pack(pady=5, anchor="w")

        botao(left, calculo_limite, "Calculate")
        botao(left, exemplo_limite, "Example")

        resultado_text_limite = ctk.CTkTextbox(right, font=font)
        resultado_text_limite.pack(fill="both", expand=True)

        img = ctk.CTkImage(Image.open("limit.png"), size=(250, 120))
        ctk.CTkLabel(right, image=img, text="").pack(pady=10)

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

        img = ctk.CTkImage(Image.open("raiz.png"), size=(250, 120))
        ctk.CTkLabel(right, image=img, text="").pack(pady=10)

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
        ctk.CTkCheckBox(left, text="Piecewise Function (separate with ';')", variable=is_piecewise_var).pack(pady=5, anchor="w")
       

        show_points_var = ctk.BooleanVar(value=False, master=left)
        ctk.CTkCheckBox(left, text="Show critical and inflection points", variable=show_points_var).pack(pady=5, anchor="w")

        botao(left, plot_grafico, "Plot")

        
        ctk.CTkButton(left, text="Import points file", command=carregar_arquivo_pontos).pack(pady=10, anchor="w")
        interpolar_var = ctk.BooleanVar(value=False,master=left)
        check_interpolar = ctk.CTkCheckBox(left, text="Interpolate curve", variable=interpolar_var)
        botao_plot_dados = ctk.CTkButton(left, text="Plot imported data", command=plotar_dados_importados)

        resultado_text_grafico = ctk.CTkTextbox(right, font=font, height=150)
        resultado_text_grafico.pack(fill="x", pady=(0, 10))

       
        frame_grafico_container = ctk.CTkFrame(right)
        frame_grafico_container.pack(fill="both", expand=True)

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

        ctk.CTkButton(left, text="When to use L'Hôpital?", command=abrir_explicacao_lhopital).pack(pady=5, anchor="w")

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

        ctk.CTkButton(left, text="What is an Integral?", command=abrir_explicacao_integral).pack(pady=5, anchor="w")

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

        img = ctk.CTkImage(Image.open("integral.png"), size=(300, 180))
        ctk.CTkLabel(right, image=img, text="").pack(pady=10)

    # ====================== ABA MANUAL =========================
    def aba_manual(self, frame):
        ctk.CTkButton(
            frame,
            text="Open DDX Manual",
            command=lambda: webbrowser.open('https://drive.google.com/file/d/1XhUZMxmc4bNwYTh5FOoXdG_eztrMHVSs/view?usp=sharing'),
            width=300
        ).pack(pady=20)


# ====================== EXECUÇÃO =========================
if __name__ == "__main__":
    initial_page = InitialPage()
    initial_page.mainloop()
