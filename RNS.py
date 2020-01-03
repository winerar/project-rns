"""
Модуль для работы с системой остаточных классов.
"""
from math import ceil, factorial

def __init__(self, modules, redundant_digits=0):
    """
    принимает modules - список модулей (list) и 
    redundant_digits - количество избыточных модулей (int)
    self.modules - система оснований (list)
    self.used_modules - используемые основания (без избыточных) (list)
    self.range - полный диапазон системы (int)
    self.used_range - используемый диапазон (int)
    """
    self.modules = list(modules)
    self.used_modules = self.modules[:len(modules) - redundant_digits]
    self.range = self.get_range(self.modules)
    self.used_range = self.get_range(self.used_modules)

def modules_mutually_simple(self):
    """
    проверяет модули на взаимную простоту
    возвращает True, если все модули взаимно просты (bool)
    возвращает False, если есть не взаимно простые модули (bool)
    """
    for i in range(len(self.modules)):
        for j in range(i + 1, len(self.modules)):
            module = self.modules[i]
            second_module = self.modules[j]
            while module != second_module:
                if second_module > module:
                    second_module -= module
                else:
                    module -= second_module
            if module != 1:
                return False
    return True

def get_range(self, modules):
    """
    вычисляет диапазон системы остаточных классов
    принимает modules - список оснований (list)
    возвращает P - произведение всех оснований (int)
    """
    P = int(1)
    for module in modules:
        P *= module
    return P

def convert_to_rns(self, dec_number):
    """
    переводит десятичное целое число в СОК с заданными основаниями
    принимает dec_number - десятичное целое число (int)
    возвращает rns_number - число в СОК (list)
    """
    rns_number = []
    for module in self.modules:
        rns_number.append(dec_number % module)
    return rns_number

def convert_to_decimal(self, rns_number, modules=None):
    """
    переводит число из СОК в десятичную систему счисления
    принимает rns_number - число в СОК (list)
    принимает modules - список модулей (list)
    возвращает целое десятичное число (int)
    """
    if modules == None:
        modules = list(self.modules)
    dec_number = int(0)
    rns_range = self.get_range(modules)
    for i in range(len(modules)):
        d = int(1)
        p = rns_range // modules[i]
        k = p % modules[i]
        while d % k != 0:
            d += modules[i]
        orthogonal_basis = int(d // k * p)
        dec_number += rns_number[i] * orthogonal_basis
    return dec_number % rns_range

def correct_single_error(self, rns_number):
    """
    обнаруживает и исправляет однократную ошибку
    переводит число из СОК в десятичную систему счисления
    принимает rns_number - число в СОК (list)
    возвращает целое десятичное число (int)
    """
    if (self.convert_to_decimal(rns_number) > self.used_range):
        for i in range(len(self.modules)):
            modules = self.modules[:i] + self.modules[i+1:]
            number = rns_number[:i] + rns_number[i+1:]
            if self.convert_to_decimal(number, modules) < self.used_range:
                return self.convert_to_decimal(number, modules)
    return self.convert_to_decimal(rns_number)

def correct_2_errors_pm(self, rns_number):
    """
    метод проекций для двух ошибок
    """
    iter = 0
    
    # проверка на наличие ошибок
    dec_number = self.convert_to_decimal(rns_number)
    if dec_number < self.used_range:
        return dec_number
    
    for i_1 in range(len(self.modules)):
        for i_2 in range(i_1 + 1, len(self.modules)):
            
            iter += 1
            # удаление i_1-го и i_2-го разрядов и вычисление проекции X
            modules = self.modules[:i_1] + self.modules[i_1+1:i_2] + self.modules[i_2+1:]
            residues = rns_number[:i_1] + rns_number[i_1+1:i_2] + rns_number[i_2+1:]
            dec_number = self.convert_to_decimal(residues, modules)
            # если проекция входит в допустимый диапазон, возвращаем её значение
            if dec_number < self.used_range:
                
                #print("Итераций: ", iter)
                
                return dec_number
            
    #print("Итераций: ", iter)
    
    return -1

def correct_2_errors_mpm(self, rns_number):
    """
    модифицированный метод проекций для двух ошибок
    """
    iter = 0
    t = 2
    
    # проверка на наличие ошибок
    Y = self.convert_to_decimal(rns_number)#, self.used_modules)
    if Y < self.used_range:
        return Y
    
    n = len(self.modules) # общее количество модулей
    o = n - len(self.used_modules)
    p = factorial(n) // factorial(n - t) // factorial(t) # количество комбинаций позиций ошибок
    f = ceil(p / (factorial(o) // factorial(o - t) // factorial(t)))
    
    for c in range(1, f+1):
        # удаление разрядов, не используемых в группе
        modules = self.modules[:(f + 1 - c) * t - t] + self.modules[(f + 1 - c) * t:]
        residues = rns_number[:(f + 1 - c) * t - t] + rns_number[(f + 1 - c) * t:]
        
        for i_1 in range(len(modules)):
            for i_2 in range(i_1 + 1, len(modules)):
                
                iter += 1
                # удаление i_1-го и i_2-го разрядов и вычисление проекции X
                mod = modules[:i_1] + modules[i_1+1:i_2] + modules[i_2+1:]
                res = residues[:i_1] + residues[i_1+1:i_2] + residues[i_2+1:]
                X = self.convert_to_decimal(res, mod)
                
                if X < self.used_range:
                    #вычисление вектора v
                    v = self.convert_to_rns(X)
                    # если расстояние меньше t = 2, возвращаем значение проекции
                    if self.hamming_distance(rns_number, v) <= t:
                        
                        #print("Итераций: ", iter)
                        
                        return X
                    else:
                        
                        break
                else:
                    
                    break
                    
    print("Итераций: ", iter)
    
    return -1

def hamming_distance(self, list_1, list_2):
    """
    возвращает расстояние Хэмминга между списками list
    возвращает -1, если списки разной длины
    """
    if len(list_1) != len(list_2):
        return -1
    distance = 0
    for i in range (len(list_1)):
        if abs(list_1[i]-list_2[i]) > 0:
            distance += 1
    return distance

def add(self, *numbers):
    s = []
    for i in range(len(self.modules)):
        si = 0
        for rns_number in numbers:
            si += rns_number[i]
        s.append(si % self.modules[i])
    return s

def multiply(self, *numbers):
    m = []
    for i in range(len(self.modules)):
        mi = 1
        for rns_number in numbers:
            mi *= rns_number[i]
        m.append(mi % self.modules[i])
    return m
