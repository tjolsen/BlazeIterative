To better improve the convergence rate, we then introduce precondition conjugate gradient in the following.

#### Precondition

The preconditioner **M** of matrix **A** is chosed such that matrix ![image](https://user-images.githubusercontent.com/29106484/61173498-d1150980-a559-11e9-9cb2-a865fdb0707b.png) has smaller condition number than matrix **A**. That is, the introducing of preconditioner can lead to the decrease of condition number, which finally results in the increase of the convergence rate. Then a question is raised, what is condition number?

#### Condition Number

Condition number is used to measure how the output changes when the input has a small change. For example, for linear equation 
**Ax = b**, if **b** has small change![image](https://user-images.githubusercontent.com/29106484/61173633-90b68b00-a55b-11e9-930a-547f7252d2e7.png), then **x** would have the change ![image](https://user-images.githubusercontent.com/29106484/61173648-afb51d00-a55b-11e9-8cf2-3bccd249908f.png), where ![image](https://user-images.githubusercontent.com/29106484/61173656-cd828200-a55b-11e9-8eb3-55f3e59f4300.png) are the error vectors.

Substitute ![image](https://user-images.githubusercontent.com/29106484/61173656-cd828200-a55b-11e9-8eb3-55f3e59f4300.png) into the linear equation, we have 

![image](https://user-images.githubusercontent.com/29106484/61173757-5bab3800-a55d-11e9-97a6-db6401919878.png)

Then our concern is: how a small change of ![image](https://user-images.githubusercontent.com/29106484/61173775-91502100-a55d-11e9-852f-d3a94e1ea4c5.png) would result in the change of ![image](https://user-images.githubusercontent.com/29106484/61173778-a0cf6a00-a55d-11e9-9f3d-2d07dbb8d8e3.png). Definitely, the smaller ![image](https://user-images.githubusercontent.com/29106484/61173778-a0cf6a00-a55d-11e9-9f3d-2d07dbb8d8e3.png) is better. 

Instead of looking the ![image](https://user-images.githubusercontent.com/29106484/61173775-91502100-a55d-11e9-852f-d3a94e1ea4c5.png) and ![image](https://user-images.githubusercontent.com/29106484/61173778-a0cf6a00-a55d-11e9-9f3d-2d07dbb8d8e3.png), let us look at the normalized version, that is, how a small change of ![image](https://user-images.githubusercontent.com/29106484/61173815-1d624880-a55e-11e9-93f5-f054890d15df.png) would lead to the change of ![image](https://user-images.githubusercontent.com/29106484/61173822-30751880-a55e-11e9-97d5-2d356825266c.png). Similarlly, it is better when ![image](https://user-images.githubusercontent.com/29106484/61173873-a5e0e900-a55e-11e9-9c2b-28b3ad11a3bf.png)
 is smaller. Then our goal is to find the upper bound of ![image](https://user-images.githubusercontent.com/29106484/61173873-a5e0e900-a55e-11e9-9c2b-28b3ad11a3bf.png).
 
Suppose ![image](https://user-images.githubusercontent.com/29106484/61173930-518a3900-a55f-11e9-8c3e-10c41c44eee5.png) is the eigenvalue of matrix **A**. Then we have ![image](https://user-images.githubusercontent.com/29106484/61173965-b3e33980-a55f-11e9-9c09-0adc0710b758.png), which implies that ![image](https://user-images.githubusercontent.com/29106484/61173972-cfe6db00-a55f-11e9-984b-e628eb1140b0.png). Similarly, we have ![image](https://user-images.githubusercontent.com/29106484/61173983-09b7e180-a560-11e9-8fb6-53017f0a56b1.png). Therefore, ![image](https://user-images.githubusercontent.com/29106484/61173992-2e13be00-a560-11e9-8f12-f5e56430463d.png).

Finally, the definition of condition number of **A** is: 

![image](https://user-images.githubusercontent.com/29106484/61178459-f97d2200-a5b2-11e9-950f-d002177750f4.png),

when matrix **A** is normal. **A** is said to be well-conditioned if the condition number is small. Furthermore, the definition of a general case is:

![image](https://user-images.githubusercontent.com/29106484/61178473-16195a00-a5b3-11e9-81b9-3d431bbf8b26.png).

#### Algorithm
After knowing the benefit of using preconditioner, let us focus on the algorithm part, which is from Wiki:

<img width="380" src="https://user-images.githubusercontent.com/29106484/61248052-f0f22c00-a717-11e9-94e8-3b0659595150.png">

#### Choices of Preconditioner
Then the last step is to choose a proper preconditioner. First we decompose matrix **A** as ![image](https://user-images.githubusercontent.com/29106484/61174340-56052080-a564-11e9-9647-48c9b51d609b.png), where **L** and **D** are strictly lower matrix and diagonal matrix, respectively. Next, we introduce some popular preconditioners in the following.

1. Jacobi preconditioning: 

**M = D**

2. Gauss-Seidel precondition: 

**M = L + D**

3. Successive over-relaxation (SOR) precondition: 

![image](https://user-images.githubusercontent.com/29106484/61174418-5a7e0900-a565-11e9-9438-7fe5894320e1.png), 

where ![image](https://user-images.githubusercontent.com/29106484/61174446-ae88ed80-a565-11e9-9639-aa646769f6d7.png) is the relaxation parameter.

4. Symmetric SOR preconditioning (SSOR): 

![image](https://user-images.githubusercontent.com/29106484/61174505-c745d300-a566-11e9-95dc-28a1888ed1f5.png)

5. Incomplete Cholesky factorization: 

![image](https://user-images.githubusercontent.com/29106484/61175824-8c9a6580-a57b-11e9-8d1c-03cb5904a9d6.png),

where matrix **K** can be computed as follow: first, we decompose matrix **A** as ![image](https://user-images.githubusercontent.com/29106484/61175770-97082f80-a57a-11e9-8d70-e8d0093ce92c.png) using Cholesky factorization, where ![image](https://user-images.githubusercontent.com/29106484/61175795-0e3dc380-a57b-11e9-8122-3b7ce6580d30.png) here is a lower triangular matrix. Next, we will compute a sparse lower triangular matrix **K**, which is close to **L**. Specifically, the elements, which are zeros in matrix **A**, are set to zeros in matrix **K**. 
