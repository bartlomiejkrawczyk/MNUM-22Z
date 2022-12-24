# Metody Numeryczne - Projekt III
```s
student: Bartłomiej Krawczyk
indeks: 310774
```
# Zadanie 1

## Treść

Ruch punktu na płaszczyźnie ($x_1$, $x_2$) jest opisany równaniami:

$$
\frac{dx_1}{dt} = x_2 + x_1 (0.3 - x_1^2 - x_2^2)
$$

$$
\frac{dx_1}{dt} = - x_1 + x_2  (0.3 - x_1^2 - x_2^2)
$$

Należy obliczyć przebieg trajektorii ruchu tego punktu w przedziale [0, 20] dla warunków początkowych: $x_1(0) = 0$, $x_2(0) = 13$.

Rozwiązanie proszę znaleźć korzystając z zaimplementowanej przez siebie w języku Matlaba w formie funkcji (możliwie uniwersalnej, czyli solwera, o odpowiednich parametrach wejścia i wyjścia) metody Dormand-Prince'a czwartego rzędu przy zmiennym kroku z szacowaniem błędu techniką pary metod włożonych (DorPri45).

## Opis algorytmu

## Program

## Porównanie rozwiązania z wykorzystaniem solwera ode45 Matlaba

<!-- Wykresy -->
- trajektoria $x_1(t)$ otrzymana dwiema metodami

- trajektoria $x_2(t)$ otrzymana dwiema metodami

- trajektoria na płaszczyźnie ($x_1$, $x_2$) (inaczej: w przestrzeni fazowej) otrzymana dwiema metodami

## Komentarz

- minimalny krok $h_{min}$

- dokładność względna

- dokładność bezwzględna

## Wykresy

- zależność długości kroku od czasu

- zależność estymaty błędu od czasu

## Wnioski

- ocena poprawności wyników

- efektywności algorytmów

- która metoda lepsza?

- która metoda szybsza?

