#Библиотеки
import math
import numpy as np
import matplotlib.pyplot as plt

# Определение констант
L = 1 #Расстояние между зеркалами
alpha_1 = np.pi/6 #Угол между нормалью к поверхности зеркала и осью Z
c = 3e8 #Скорость света
dT_0_squared = 1 #d(|T0|^2)/df
T_0_squared = 0.5 #|T0|^2
omega = 2 * np.pi #Циклическая частота

#Все параметры из первой серии линий
#Параметры для эпсилон1
f_epsi1 = np.array([132.366, 132.433]) #Частота
f_width1 = 0.05 * 10**3 #Ширина df
epsi01 = 10**6 #Epsi1

#Параметры для эпсилон3
f_epsi3 = np.arange(132.3, 132.45, 0.05) #Частота для epsi3: [132.3, 132.35, 132.4, 132.45]
f_width3 = np.arange(0.03 * 10**6, 0.06 * 10**6, 0.01 * 10**6) #Ширина df: [10, 8, 6, 4] * 10^6
epsi03 = np.arange(10 * 10**6, 3* 10**6, -2* 10**6) #Epsi3

# Определение функций
#Лоренцевы контуры
#Первая цифра - какой эпсилон, вторая цифра - количество штрихов
def get_epsi_11(f):
	arr = []
	for i in range(len(f)):
		value = 0
		for j in range(len(f_epsi1)):
			value += (epsi01 * (f_epsi1[j]**2 - f[i]**2) * f_epsi1[j] * f_width1) / ((f_epsi1[j]**2 - f[i]**2)**2 + f_width1**2 * f[i]**2)
		arr.append(1 + value)
	return arr

def get_epsi_12(f):
	arr = []
	for i in range(len(f)):
		value = 0
		for j in range(len(f_epsi1)):
			value += (epsi01 * f_width1**2 * f_epsi1[j] * f[i]) / ((f_epsi1[j]**2 - f[i]**2)**2 + f_width1**2 * f[i]**2)
		arr.append(1 + value)
	return arr

def get_epsi_31(f):
	arr = []
	for i in range(len(f)):
		value = 0
		for j in range(len(f_epsi3)):
			value += (epsi03[j] * (f_epsi3[j]**2 - f[i]**2) * f_epsi3[j] * f_width3[j]) / ((f_epsi3[j]**2 - f[i]**2)**2 + f_width3[j]**2 * f[i]**2)
		arr.append(1 + value)
	return arr

def get_epsi_32(f):
	arr = []
	for i in range(len(f)):
		value = 0
		for j in range(len(f_epsi3)):
			value += (epsi03[j] * f_width3[j]**2 * f_epsi3[j] * f[i]) / ((f_epsi3[j]**2 - f[i]**2)**2 + f_width3[j]**2 * f[i]**2)
		arr.append(1 + value)
	return arr

def get_kz32(epsi31, epsi32):
	arr = []
	for i in range(len(epsi32)):
		arr.append(omega * epsi32[i] / c * 2 * np.sqrt(epsi31[i] * np.cos(alpha_1)))
	return arr

def get_alpha_3(epsi31, epsi32):
	arr = []
	for i in range(len(epsi31)):
		arr.append(math.asin((epsi31[i] * math.sin(math.radians(alpha_1))) / epsi32[i]))
	return arr

def T_sample_squared(f, f0):
	arr = []
	for i in range(len(f)):
		arr.append((T_0_squared + ((np.sqrt(epsi3_arr[i] * np.cos(alpha_1) / epsi1_arr[i] * np.cos(alpha_3[i]))) - 1) * dT_0_squared * (f[i] - f0[i])) * np.exp(-2 * kz32[i] * L))
	return arr

#Остальные параметры
width = 0.06 * 10**-3 #Ширина для центральной частоты
f_values_center = np.linspace(132.3 - 3 * width, 132.3 + 3 * width, 60) #Центральная частота линии: (f0 - 3df), (f0 + 3df)
f_values = np.array([f_values_center[i - 1] - 3 * width + i * width / 10 for i in range(1, 61)]) #Массив частот для как переменная для основной формулы и эпсилон

epsi1_arr = get_epsi_11(f_values) #Эпсилон1 штрих
epsi12_arr = get_epsi_12(f_values) #Эпсилон1 с двумя штрихами
epsi3_arr = get_epsi_31(f_values) #Эпсилон3 с одним штрихом
epsi32_arr = get_epsi_32(f_values) #Эпсилон3 с двумя штрихами
kz32 = get_kz32(epsi3_arr, epsi32_arr) #Волновой вектор
alpha_3 = get_alpha_3(epsi3_arr, epsi32_arr) #Альфа 3 через акрсинус
T_sample_squared_values = T_sample_squared(f_values, f_values_center) #Коэффициент пропускания

# Построение графика
plt.plot(f_values, T_sample_squared_values)
plt.xlabel('f')
plt.ylabel('|T_sample|^2')
plt.title('Зависимость |T_sample|^2 от f')
plt.grid(True)
plt.show()