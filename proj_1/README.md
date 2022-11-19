# Metody Numeryczne - Projekt I
```s
student: Bartłomiej Krawczyk
indeks: 310774
```


# Zadanie 1

## Treść

Napisać uniwersalną procedurę w Matlabie o odpowiednich parametrach wejścia i wyjścia (solwer), rozwiązującą układ $n$ równań liniowych $Ax = b$, wykorzystując podaną metodę.

> Nie sprawdzać w procedurze, czy dana macierz $A$ spełnia wymagania stosowalności metody.

Obliczyć błąd rozwiązania $\varepsilon = ∥A\tilde{x} − b∥_2$ (skorzystać z funkcji `norm` Matlaba).

Proszę zastosować następnie swoją procedurę w programie do rozwiązania obydwu (jeśli można) lub jednego z układów równań dla podanych niżej macierzy $A$ i wektorów $b$, przyjmując: $n = 5, 10, 25, 50, 100, 200$.

Metoda: faktoryzacji $LDL^T$

Proszę wykonać wykres (wykresy) zależności błędu $\varepsilon$ od liczby równań $n$.

## Rozwiązanie



$$
\begin{align}
y = a
\end{align}
$$

```matlab
function x = solve_linar_equation_with_triangular_matrix_upper(A, b)

    [n, ~] = size(A);
    x = zeros(n, 1);
    
    for k = n : -1 : 1
        x(k, 1) = b(k, 1);
        
        for j = k + 1 : n
            x(k, 1) = x(k, 1) - A(k, j) * x(j, 1);
        end

        x(k, 1) = x(k, 1) / A(k, k);
    end
end
```


# Zadanie 2

## Treść

Napisać uniwersalną procedurę w Matlabie o odpowiednich parametrach wejścia i wyjścia, rozwiązującą układ $n$ równań liniowych $Ax = b$, wykorzystując metodę **iteracyjną Jacobiego**.

> Nie sprawdzać w procedurze, czy dana macierz $A$ spełnia wymagania stosowalności metody.

Jej parametry wejściowe powinny zawierać m.in. wartość graniczną $\delta$ błędu między kolejnymi przybliżeniami rozwiązania, liczonego jako **norma euklidesowa z ich różnicy** (skorzystać z funkcji `norm` Matlaba). Przyjąć jako kryterium stopu warunek: δ = 10−8 ≜ 1e − 8.

$$
\delta = 10^{-8} \triangleq 1e-8
$$

Proszę zastosować tę procedurę do rozwiązania właściwego układu równań spośród przedstawionych poniżej dla $n = 5, 10, 25, 50, 100, 200$.

Proszę sprawdzić dokładność rozwiązania licząc także błąd $\varepsilon$ i dla każdego układu równań wykonać rysunek zależności tego błędu od liczby równań $n$. Jeśli był rozwiązywany ten sam układ równań, co w p. 1, proszę porównać czasy obliczeń dla różnych algorytmów i wymiarów zadań.

## Rozwiązanie

# Zadanie 3

## Treść

Dla podanych w tabeli danych pomiarowych (próbek) metodą najmniejszych kwadratów należy wyznaczyć funkcję wielomianową y = f(x) (tzn. wektor współczynników) najlepiej aproksymującą te dane.

| $x_i$ | $y_i$   |
|-------|---------|
| -10   | -42.417 |
| -8    | -23.440 |
| -6    | -11.160 |
| -4    | -4.128  |
| -2    | -0.725  |
| 0     | 0.942   |
| 2     | -2.069  |
| 4     | -3.908  |
| 6     | -4.705  |
| 8     | -5.438  |
| 10    | -3.578  |

Proszę przetestować wielomiany stopni: $3, 5, 7, 9, 10$. Kod aproksymujący powinien być uniwersalną procedurą w Matlabie o odpowiednich parametrach wejścia i wyjścia.

W sprawozdaniu proszę przedstawić na rysunku otrzymaną funkcję na tle danych (funkcję aproksymującą proszę próbkować przynajmniej 10 razy częściej niż dane).

Do rozwiązania zadania najmniejszych kwadratów proszę wykorzystać najpierw **układ równań normalnych**, a potem **rozkład SVD**.

Do rozwiązywania układu równań i dekompozycji użyć solwerów Matlaba.. Porównać efektywność obydwu podejść.

Do liczenia wartości wielomianu użyć funkcji `polyval`.

Proszę obliczyć błąd aproksymacji w dwóch normach: euklidesowej oraz maksimum (nieskończoność). W obydwu przypadkach skorzystać z funkcji `norm` Matlaba.

## Rozwiązanie

