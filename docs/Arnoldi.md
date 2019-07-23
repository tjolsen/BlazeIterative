Arnoldi iteration is an iterative method to approximately find the eigenvalues and eigenvectors of matrices by generating an orthogonal basis of Krylov subspace. 

#### Krylov Subspace
Given a matrix **A** which is m x m, and an initial vector **b**, which is m x 1, we can construct the n-th order Krylov subspace:

![image](https://user-images.githubusercontent.com/29106484/61185593-d2f5d000-a620-11e9-8dde-e29f0db78fb0.png).

The advantage of constructing the Krylov subspace is that we can compute the matrix-vector product instead of cope with matrix **A** directly, which is particularly useful when **A** is large, since the cost of matrix-vector product is relatively cheaper.

Suppose ![image](https://user-images.githubusercontent.com/29106484/61185779-b9558800-a622-11e9-9412-a97c576c129f.png), where **H** and **Q** are m x m matrices, **H** is an upper Hessenberg matrix, which has zero entries below the first subdiagonal, and  **Q** is a unitary matrix with ![image](https://user-images.githubusercontent.com/29106484/61186005-8791f080-a625-11e9-8acb-f37c0a184472.png). With the aid of Krylov subspace, we can reduce the matrix **A** to an upper Hessenberg matrix. 

First, we can rewritten the expression of matrix **A** as **AQ = QH**, shown as follows:

![image](https://user-images.githubusercontent.com/29106484/61188119-34796700-a640-11e9-901f-ddf08f36dce2.png),

where n is the dimention of Krylov subspace as mentioned. Then we write a new expression of matrix **A**, which is part of the above equation, as ![image](https://user-images.githubusercontent.com/29106484/61188186-26781600-a641-11e9-9e90-1dffd1996170.png), where ![image](https://user-images.githubusercontent.com/29106484/61188198-43ace480-a641-11e9-98f5-7e0fd32f2f09.png) and 

![image](https://user-images.githubusercontent.com/29106484/61188176-05172a00-a641-11e9-94a6-560705a0619b.png).

Calculated the n-th columns of both sides, we have ![image](https://user-images.githubusercontent.com/29106484/61188268-49ef9080-a642-11e9-95bb-2411e1d2f0f4.png). Therefore, we can express  ![image](https://user-images.githubusercontent.com/29106484/61188294-bbc7da00-a642-11e9-82b8-47fb289e55a6.png) as

![image](https://user-images.githubusercontent.com/29106484/61188293-b074ae80-a642-11e9-9879-471aae3062f9.png),

which can be regarded as a iterative process.

#### Algorithm
The algorithm is shown as follow, which is from wiki:

<img width="360" alt="arnoldi" src="https://user-images.githubusercontent.com/29106484/61188359-7a83fa00-a643-11e9-84dd-237d41a29ecf.png">.

Note that ![image](https://user-images.githubusercontent.com/29106484/61189180-8a094000-a64f-11e9-9a4d-2cb2add138a7.png), which can be computed using the above algorithm. By removing the last row of matrix ![image](https://user-images.githubusercontent.com/29106484/61189338-99898880-a651-11e9-9ac7-53162e137e59.png), we have the n x n matrix ![image](https://user-images.githubusercontent.com/29106484/61189306-2c75f300-a651-11e9-9f52-93d045929020.png), which has the same eigenvalues as matrix **A**. That is how we reduce the matrix **A** to an upper Hessenberg matrix.
