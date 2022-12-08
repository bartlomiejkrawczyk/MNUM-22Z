# Metody Numeryczne - Projekt II
```s
student: Bartłomiej Krawczyk
indeks: 310774
```
# Zadanie 1

## Treść

Proszę znaleźć wszystkie pierwiastki funkcji
$f(x) = 2.1 - 2x - e^{-x/2}$ w przedziale $[10, -10]$ używając dla każdego zera:

a) własnego solwera z implementacją metody **siecznych**

b) podanego na stronie przedmiotu solwera `newton.m` z implementacją metody Newtona

## Rozwiązanie a)

W metodzie siecznych prowadzimy sieczną zawsze między dwoma ostatnio wyznaczonymi punktami. Nie dbamy przy tym o zachowanie przedziału izolacji pierwiastka.

Wyznaczamy rozwiązanie układu równań - prostej $y=0$ oraz prostej przechodzącej przez dwa ostatnio wyznaczone punkty:
$$
\begin{cases}
y = 0 \\
(y - f(x_n))(x_{n-1} - x_n) - (f(x_{n-1}) - f(x_n))(x - x_n) = 0\\
\end{cases}
$$
gdzie rozwiązaniem $x$ będzie $x_{n+1}$. Otrzymujemy:

$$
x_{n+1} = x_n - \frac{f(x_n)(x_n - x_{n-1})}{f(x_n) - f(x_{n-1})} = \frac{x_{n-1}f(x_n) - x_n f(x_{n-1})}{f(x_n) - f(x_{n-1})}
$$

## Rozwiązanie b)


# Zadanie 2

## Treść

Używając metody **Müllera MM1** proszę znaleźć wszystkie pierwiastki wielomianu czwartego stopnia:

$$
f(x) = a_4x^4 + a_3x^3 + a_2x^2 + a_1x + a_0
$$


$$
\begin{bmatrix}
a_4 & a_3 & a_2 & a_1 & a_0\\
\end{bmatrix}
=
\begin{bmatrix}
-1 & 1.5 & 1.5 & 0.5 & 1\\
\end{bmatrix}
$$
