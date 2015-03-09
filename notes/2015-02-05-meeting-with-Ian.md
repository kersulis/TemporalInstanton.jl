# Meeting with Ian Hiskens
_February 5, 2015_

## Updates

We've both been implementing Dan Bienstock's "trust-region subproblem" procedure, as described in his paper and emails.

Non-negativity constraints: ignore temporarily. Translation will result in $ Ax \geq k $. One paper referenced by Dan made progress by assuming polytope intersection points are outside the region of interest. See first page of "Polynomial Solvability": the region of interest contains no points for which $a_i^\top x = b_i$ and $a_j^\top x = b_j$.

Does $G$ being positive semidefinite cause problems later on?

Keep all variables in formulation: renewable generation and voltage angles.

For now, consider every decision variable to be independent. In the long term, look at auto- and cross- correlation between wind generation variables. Can we capture this correlation between times and between sites?

Decide how to order variables. Need to get model set out and get all information into Dan's solution algorithm. What comes out of this is a set of points that need to be enumerated and checked. 

Scaling each variable to obtain norm constraint: Norm constraint only involves subset of variables. Could introduce new variables by multiplying by square root of Q. Use diagonal matrix with powers of $\tau$, adjust $A$ by same matrix. This is just part of conditioning the model prior to solving. If we choose a different slack node for each line, the norm constraint has one variable. But then we have a different set of nodes when considering each line.

## Steps forward

1. Get this working for **one line**.

2. Map our problem into the structure given by Dan!

3. Finish Dan's pipeline.

4. This week: translate problem structure into what Dan sent us. What is $G$ in terms of $Y$ and other network parameters?

Could build deviation into $Ax=b$, or could keep generation by itself.

When working with temperature constraint, need actual angles -- not angles corresponding to forecast!

Note: forecast is fully know a priori. We can fix base power flows at all nodes and time steps using forecast information.

If forecast is entirely accurate, we can compute what the line temperature will be. This can give us preliminary filtering. Are there any lines that are forecast to be hot? Suppose the forecast will cause a line to overheat: 

1. Structure: mapping through to Dan's formulation.

2. Solution: does the output make sense?

3. Application: how can we use this?