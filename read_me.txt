
Purpose:
--------

This package implements the "revised simplex" algorithm for linear programs in standard form.

It embraces all type of constraints, in other words, if an initial basic feasible solution is

trivial, then the algorithm finds it and then solves the LP with revised simplex.Otherwise,if

an initial basic feasible solution is not trivial due to equality or greater than or equal to

constraints, then it applies "Two  Phase" approach, each phase of the method, revised simplex

algorithm is used to find the optimal solution. 

How to use:
-----------

Due to object-oriented design, it is very easy to employ this package in your linear programs.

In order to illustrate how you can use it, '_test1.py' and '_test2.py' are created as an example

cases. 


* '_test1.py' file uses the example 7.2.1 in Metin Türkay's INDR 501 Lecture Notes. This is an 

example of how to solve your LP when an initial basic feasible solution is trivial. As a result,

the algorithm finds the optimal solution in a short amount of time.


* '_test2.py' file uses the example 6.2.1 in Metin Türkay's INDR 501 Lecture Notes. This is an 

example of how to solve your LP when an initial basic feasible solution is not trivial. As a result,

the algorithm finds the optimal solution in a short amount of time.


