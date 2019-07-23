We will introduce conjugate gradient method in the following. Before we introduce the iterative algorithms, let us look at the direct method, which enables us gain a better understanding of the iterative algorithms. 

#### Direct Method of Conjugate Gradient

Suppose that 

![conjugate vectors p](https://user-images.githubusercontent.com/29106484/61144125-a4012200-a499-11e9-94a2-1b7511091de8.png),

where 

![what are conjugate vectors](https://user-images.githubusercontent.com/29106484/61144273-140fa800-a49a-11e9-873c-0da839707a09.png),

is a set of conjugate vectors with respect to matrix A. 

Let ![x*](https://user-images.githubusercontent.com/29106484/61145178-62be4180-a49c-11e9-8fda-db6b4d99ebdd.png) denote the solution. We can express it as

![x* expression](https://user-images.githubusercontent.com/29106484/61144380-5933da00-a49a-11e9-9661-61f007c730fe.png).

Therefore, we have

![x* second expressin](https://user-images.githubusercontent.com/29106484/61144490-a021cf80-a49a-11e9-9cdc-1c7770f6ff0b.png).

Therefore, we have 

![alpha](https://user-images.githubusercontent.com/29106484/61144604-e37c3e00-a49a-11e9-8927-23b379fce2c9.png),

which implies that we can get  ![alpha_i](https://latex.codecogs.com/gif.latex?%5Calpha_%7Bi%7D)  without knowing ![x*](https://user-images.githubusercontent.com/29106484/61145178-62be4180-a49c-11e9-8fda-db6b4d99ebdd.png)
.

Substitute  ![alpha_i](https://latex.codecogs.com/gif.latex?%5Calpha_%7Bi%7D)  into ![x*](https://user-images.githubusercontent.com/29106484/61145178-62be4180-a49c-11e9-8fda-db6b4d99ebdd.png)
, we have

![x* final expression](https://user-images.githubusercontent.com/29106484/61144684-2b02ca00-a49b-11e9-9d1a-c49f910a1570.png).

Observing the above expression, we find that we do not need to calcaulate matrix inversion. Furthermore, the expression can be regarded as iterative process, wherein the n-th term is added at the n-th iteration.

#### Basic Iterative Method of Conjugate Gradient

As we mentioned above, the direct menthod is costly when n is large. To avoid such cost, we dynamic gererate the conjugate vectors instead of finding them via direct method.

Let ![x1](https://user-images.githubusercontent.com/29106484/61148338-27c00c00-a4a4-11e9-81dc-bf1ddb926b33.png) denote the initial guess, the update rule is shown as follow:

At k-th iteration:

![update rule](https://user-images.githubusercontent.com/29106484/61148860-568ab200-a4a5-11e9-9fa4-be12bbc9db0c.png)

where ![image](https://user-images.githubusercontent.com/29106484/61148953-918ce580-a4a5-11e9-801c-62e2734dbd04.png) is the redisual at the k-th iteration.

Then, after n-th iteration, we have ![image](https://user-images.githubusercontent.com/29106484/61160644-d6287900-a4c5-11e9-8247-e1b138aec760.png).

The proof is shown as follow:
First, we express ![image](https://user-images.githubusercontent.com/29106484/61149957-f0535e80-a4a7-11e9-94cd-de0cea69801f.png) as

![expression of xs](https://user-images.githubusercontent.com/29106484/61160320-a0cf5b80-a4c4-11e9-9fb6-e9351eee5be2.png).

Therefore, we have

![image](https://user-images.githubusercontent.com/29106484/61160613-ab3e2500-a4c5-11e9-8b7c-7ce87a6df824.png)


Therefore, we have

![image](https://user-images.githubusercontent.com/29106484/61160738-46cf9580-a4c6-11e9-9a64-3453c37abe20.png)

and we can express ![alpha_k](https://user-images.githubusercontent.com/29106484/61150457-188f8d00-a4a9-11e9-9456-f0b13b743838.png)
 as

![expression of alpha_k](https://user-images.githubusercontent.com/29106484/61150418-044b9000-a4a9-11e9-9602-aac99c67a0f5.png)

Also, we have

![image](https://user-images.githubusercontent.com/29106484/61151653-0c58ff00-a4ac-11e9-9df4-da8d60a4efe4.png)

Therefore, we have

![image](https://user-images.githubusercontent.com/29106484/61151743-44f8d880-a4ac-11e9-9431-441bf473ae31.png)

Finally, we can express ![alpha_k](https://user-images.githubusercontent.com/29106484/61150457-188f8d00-a4a9-11e9-9456-f0b13b743838.png) as

![image](https://user-images.githubusercontent.com/29106484/61151918-ac168d00-a4ac-11e9-9ba8-059c702278b9.png).


#### Improved Iterative Method of Conjugate Gradient

However, the above basic iterative method is still computationally expensive due to that it has to store all previous redisual vectors. A promising approach to avoid such cost is to generate a new conjugate vector by only using the previous one. Specifically, we determine the new conjugate vector by the following formular:

![image](https://user-images.githubusercontent.com/29106484/61161450-569ca900-a4c9-11e9-9778-8752fdabb34c.png)

where 

![image](https://user-images.githubusercontent.com/29106484/61161826-52718b00-a4cb-11e9-86a4-df2d4de6b4e2.png)

and

![beta_k](https://user-images.githubusercontent.com/29106484/61161555-e4789400-a4c9-11e9-8897-ff075b0da674.png).

Note that the proof of the above ![image](https://user-images.githubusercontent.com/29106484/61163774-f4967080-a4d5-11e9-9651-664ebe6ebe6e.png) expression is shown as follow:

Since ![image](https://user-images.githubusercontent.com/29106484/61163742-c7e25900-a4d5-11e9-871e-efc6479af677.png),

we have 

![image](https://user-images.githubusercontent.com/29106484/61163547-4d650980-a4d4-11e9-9b10-f5237d39cc42.png). 

Therefore, the numertor and denominator can be expressed as

![image](https://user-images.githubusercontent.com/29106484/61163638-e1cf6c00-a4d4-11e9-9aa3-38b640f2247e.png)

and

![image](https://user-images.githubusercontent.com/29106484/61163704-7b971900-a4d5-11e9-9a74-2ddb7e88806b.png),

respectively. Therefore, we can express ![image](https://user-images.githubusercontent.com/29106484/61163774-f4967080-a4d5-11e9-9651-664ebe6ebe6e.png) as the above.

Finally, the algorithm is shown as follow:
Let ![x1](https://user-images.githubusercontent.com/29106484/61148338-27c00c00-a4a4-11e9-81dc-bf1ddb926b33.png) denote the initial guess, ![image](https://user-images.githubusercontent.com/29106484/61161991-125ed800-a4cc-11e9-8b6a-7e0485e02552.png).

At k-th iteration:

![image](https://user-images.githubusercontent.com/29106484/61162762-d75ea380-a4cf-11e9-8f35-5131b75ce338.png)
