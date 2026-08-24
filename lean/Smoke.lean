import Mathlib.Tactic.Ring
import Mathlib.Data.Int.GCD

-- prueba de humo: que el entorno rapido compila y que `decide` cierra
-- aritmetica modular sobre enteros concretos, que es todo lo que los
-- seis ladrillos necesitan.
example : (-27 : ℤ) - 2 * (-18) = 9 := by decide
example : (9 : ℤ) % 6 = 3 := by decide
example (a b : ℤ) : (a + b) - 2 * (a + b) = -(a + b) := by ring
