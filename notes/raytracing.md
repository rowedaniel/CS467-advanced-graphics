
# numerical derivatives:
  - just use limit definition

## 2D: 
df   f(xbar) - f(x)    f(x+h)-f(x)
-- = -------------- = --------------
dx      xbar - x            h

## 3D:
df   f(xbar, y) - f(x,y)
-- = --------------
dx      xbar - x

df   f(x, ybar) - f(x,y) 
-- = --------------
dx      ybar - y

# normal vectors:
## 2D:
<$$frac{df}{dx}$$, $$frac{df}{dy}$$>
[df/dx, df/dy]
## 3D
<$$frac{df}{dx}$$, $$frac{df}{dy}$$, $$frac{df}{dz}$$>
[df/dx, df/dy, df/dz]


# problems for transforming from object space to world space

Using the above method, it's easy to get the normal vector in object space, by calculating df/dx, df/dy, df/dz (either
numerically or otherwise.)

It's harder to get it for object space. You can't just multiply by the view matrix: asymmetric scaling in the
transformation will distort the vector.

In order to get the normal in world space, you have to chain rule: (p = partialj)
dF   pf   dx   pf   dy
-- = -- * -- + -- * --
du   px   du   py   du

(and likewise for dv)

This ends up being the following matrix multiplication:

| dx/du  dy/du | (df/dx)
| dx/dv  dy/dv | (df/dy)

(The matrix is the Jacobian, and constant per equation. Calculate it once, and then multiply)
