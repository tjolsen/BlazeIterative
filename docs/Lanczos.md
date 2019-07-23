Different from Arnoldi, here we suppose matrix **A** is an n x n Hermitian matrix and we plan to find the m most valueable eigenvalues. In the following, we will introduce Lanczos iterative algorithm.

Let matrix **V** be an n x m matrix, which can be written as ![image](https://user-images.githubusercontent.com/29106484/61256376-3ae80b80-a732-11e9-9983-13fed202c8f7.png), where ![image](https://user-images.githubusercontent.com/29106484/61256466-81d60100-a732-11e9-82a7-e7e1b5862d95.png) are orthonormal Lanczos vectors. Let matrix **T** be a m x m tridiagonal real symmetric matrix with ![image](https://user-images.githubusercontent.com/29106484/61256853-037a5e80-a734-11e9-83c5-4d93c08ffaf6.png), which can be also written as

![image](https://user-images.githubusercontent.com/29106484/61256777-a1b9f480-a733-11e9-8f14-2ada201c7b02.png).

Note that ![image](https://user-images.githubusercontent.com/29106484/61256853-037a5e80-a734-11e9-83c5-4d93c08ffaf6.png) can be expressed as **AV = VT**. Comparing the j-th column of both sides, we have

![image](https://user-images.githubusercontent.com/29106484/61257040-e72af180-a734-11e9-8efb-91602ca90802.png). 

Furthermore, we have 

![image](https://user-images.githubusercontent.com/29106484/61257078-0e81be80-a735-11e9-92ee-bf14df8a64ab.png)

since ![image](https://user-images.githubusercontent.com/29106484/61257340-40475500-a736-11e9-97f2-a7772a4f9bb4.png).

In addition, we have

![image](https://user-images.githubusercontent.com/29106484/61257170-7df7ae00-a735-11e9-89cb-05cecb6fe48d.png) 

by multiplying ![image](https://user-images.githubusercontent.com/29106484/61257126-44bf3e00-a735-11e9-8fdb-5328125682d0.png) at both sides of ![image](https://user-images.githubusercontent.com/29106484/61257040-e72af180-a734-11e9-8efb-91602ca90802.png) since ![image](https://user-images.githubusercontent.com/29106484/61257206-a54e7b00-a735-11e9-905f-4b4b77a9d429.png), ![image](https://user-images.githubusercontent.com/29106484/61257210-aed7e300-a735-11e9-82d4-ba2b2dbb696a.png), and ![image](https://user-images.githubusercontent.com/29106484/61257217-b8614b00-a735-11e9-91a7-ac0053d59c6b.png) are orthogonal vectors.

However, in order to guarantee numerical stability, ![image](https://user-images.githubusercontent.com/29106484/61257446-a3d18280-a736-11e9-926c-20ab5c5383d1.png) and ![image](https://user-images.githubusercontent.com/29106484/61257458-adf38100-a736-11e9-9d50-2c31da4f44fb.png) should add additional terms. Therefore, they are expressed as follows:

1. ![image](https://user-images.githubusercontent.com/29106484/61257515-07f44680-a737-11e9-8791-19c3e4f04f3f.png),

since ![image](https://user-images.githubusercontent.com/29106484/61257210-aed7e300-a735-11e9-82d4-ba2b2dbb696a.png) and ![image](https://user-images.githubusercontent.com/29106484/61257217-b8614b00-a735-11e9-91a7-ac0053d59c6b.png) are orthogonal vectors.

2. ![image](https://user-images.githubusercontent.com/29106484/61257590-573a7700-a737-11e9-8e94-4389d4b0fde1.png)

since ![image](https://user-images.githubusercontent.com/29106484/61257206-a54e7b00-a735-11e9-905f-4b4b77a9d429.png) has unit norm.

The algorithm is shown as follow, which is from wiki:

<img width="792" alt="Lanczos" src="https://user-images.githubusercontent.com/29106484/61257752-fb242280-a737-11e9-9297-d238aabadb3b.png">

