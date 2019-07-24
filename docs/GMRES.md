Previously, we introduce several iterative algorims to solve the linear equation **Ax = b**,
where **A** is n x n symmetric matrix. What if **A** is non-symmetric matrix? To address the issue, we apply the generalized minimal residual method (GMRES), which will be introduced in the following.

The basic idea of GMRES is to construct the approximate solution ![x_k](https://user-images.githubusercontent.com/29106484/61744118-245a3980-ad5c-11e9-80c7-157a8c3daca1.png), where ![image](https://user-images.githubusercontent.com/29106484/61745044-4d7bc980-ad5e-11e9-82a0-b4004d4917aa.png). Note that ![image](https://user-images.githubusercontent.com/29106484/61744521-21ac1400-ad5d-11e9-8f58-eb549be1cb4c.png) is the k-dimension Krylov subspace,  ![image](https://user-images.githubusercontent.com/29106484/61744592-4bfdd180-ad5d-11e9-97d0-730d298e1dcf.png) is the initial guess, and ![image](https://user-images.githubusercontent.com/29106484/61744636-69cb3680-ad5d-11e9-9231-f768d5ccfb82.png) is the initial residual.

#### Change Target Function

Recall that our goal is to find the solution that minimizes the residual, shown as follow: ![image](https://user-images.githubusercontent.com/29106484/61746931-68e8d380-ad62-11e9-9f2d-3213051963e8.png)

Let ![image](https://user-images.githubusercontent.com/29106484/61745398-1eb22300-ad5f-11e9-81b2-02cc4477e0df.png), where ![image](https://user-images.githubusercontent.com/29106484/61745672-aac44a80-ad5f-11e9-902c-f665e9e98f52.png) is k x k unitary matrix and can be computed by using Arnoldi iteration.

Therefore, we can rewrite the function to be minimized as:

![image](https://user-images.githubusercontent.com/29106484/61752469-65a91400-ad71-11e9-883c-a3931e63faaa.png),

where (a) equation holds since 

![image](https://user-images.githubusercontent.com/29106484/61752264-abb1a800-ad70-11e9-8be7-8e9944995aff.png),

(b) equation holds according to the Arnoldi equation ![image](https://user-images.githubusercontent.com/29106484/61751536-5a081e00-ad6e-11e9-8492-547903420c18.png),

(c) equation holds since ![image](https://user-images.githubusercontent.com/29106484/61752522-8d987780-ad71-11e9-8ada-f21af4583171.png) is the first column of (k+1) x (k+1) identity matrix,

and (d) equation holds since ![image](https://user-images.githubusercontent.com/29106484/61752567-bcaee900-ad71-11e9-8b11-cb14e9c2c6e5.png) is orthonormal.

#### Algorithm
In summary, the algorithm is shown as follow

At k-th iteration:
1. Compute ![image](https://user-images.githubusercontent.com/29106484/61754329-51b4e080-ad78-11e9-9865-c0441b5f8f46.png) using Arnoldi iteration;
2. Find ![image](https://user-images.githubusercontent.com/29106484/61754386-9771a900-ad78-11e9-83af-387414cece3f.png) that minimizes ![image](https://user-images.githubusercontent.com/29106484/61754736-03084600-ad7a-11e9-99d2-b8246d3c4999.png);
3. Compute ![image](https://user-images.githubusercontent.com/29106484/61754441-d1db4600-ad78-11e9-84d4-e445f8f777f7.png);
4. Repeat if residual is not small yet.

#### Solve Least Square Problem
In step 3 of above algorithm, we plan to find the ![image](https://user-images.githubusercontent.com/29106484/61754386-9771a900-ad78-11e9-83af-387414cece3f.png) that minimizes ![image](https://user-images.githubusercontent.com/29106484/61754736-03084600-ad7a-11e9-99d2-b8246d3c4999.png). In the following, let us take a close look at the detailed implemention. 

To solve the least square problem 

![image](https://user-images.githubusercontent.com/29106484/61755151-ebca5800-ad7b-11e9-85e6-a51efdcb233f.png), 

we adopt QR decomposition, shwon as follow

![image](https://user-images.githubusercontent.com/29106484/61755562-d3f3d380-ad7d-11e9-8bfe-02a52d9e0c8f.png),

where ![image](https://user-images.githubusercontent.com/29106484/61755585-f38afc00-ad7d-11e9-9b59-d9a5075f39ce.png) is (k+1) x (k+1) orthogonal matrix and ![image](https://user-images.githubusercontent.com/29106484/61755618-0ef60700-ad7e-11e9-8329-b6f8522df30d.png) is (k+1) x k upper triangular matrix.

The difficulty is that we expect to update the decomposition of ![image](https://user-images.githubusercontent.com/29106484/61755723-77dd7f00-ad7e-11e9-868e-e10b7d45aca7.png) cheaply at each step of Arnoldi iteration. That is,

![image](https://user-images.githubusercontent.com/29106484/61755894-33061800-ad7f-11e9-894b-e061729a24ae.png),

where ![image](https://user-images.githubusercontent.com/29106484/61755943-6cd71e80-ad7f-11e9-946f-e6cb7bbf13a4.png).

To proceed, we apply Given rotations.

