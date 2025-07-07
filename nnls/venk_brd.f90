!> Butler–Reeds–Dawson regularization translated from Python.
!> Arrays are expected to be in column-major (Fortran) order.
!> The algorithm iteratively solves for c⃗ by Newton's method and
!> updates α according to the BRD prescription.
      subroutine venk_brd(k0,m,n,mvec,fvec,alpha_out,initial_alpha,tol,maxiter)
      implicit none
      integer,intent(in) :: m,n,maxiter
      double precision,intent(in) :: k0(m,n),mvec(m),initial_alpha,tol
      double precision,intent(out) :: fvec(n),alpha_out

      double precision :: alpha,alpha_new,sqrt_n
      double precision :: norm_grad
      integer :: i,j,iter
      integer :: info
      double precision,dimension(m) :: c,grad,delta_c
      double precision,dimension(m,m) :: gmat,hmat
      double precision,dimension(n) :: ktc
      double precision,dimension(m,m) :: hcopy
      integer,dimension(m) :: ipiv

      external dgesv

      sqrt_n = sqrt(dble(m))
      alpha = initial_alpha
      c = 1d0 ! start Newton at all-ones vector

      do iter=1,maxiter
         ! 01 g(c) = sum_k k(:,k) k(:,k)^T for positive (K^T c)_k
         call compute_gmat(k0,m,n,c,gmat)
         ! 02 ∇χ(c) = g(c)·c + α c − m
         call compute_grad(gmat,c,mvec,alpha,grad)
         norm_grad = sqrt(sum(grad**2))
         if (norm_grad/sqrt(sum(mvec**2)) < 1d-8) exit
         ! 03 H(c) = diag(g(c)) + α I
         call compute_hessian(gmat,alpha,hmat)
         ! 04 Newton step using LAPACK dgesv
         hcopy = hmat
         delta_c = grad
         call dgesv(m,1,hcopy,m,ipiv,delta_c,m,info)
         c = c - delta_c
         ! check convergence of gradient
         call compute_grad(gmat,c,mvec,alpha,grad)
         if (sqrt(sum(grad**2))/sqrt(sum(mvec**2)) < 1d-8) exit
         ! 06 α update
         alpha_new = sqrt_n / sqrt(sum(c**2))
         if (abs(alpha_new-alpha)/alpha < tol) exit
         alpha = alpha_new
      end do

      ! 05 recover f = max(0, K₀ᵀ·c)
      ktc = matmul(transpose(k0), c)
      do i=1,n
         if (ktc(i) > 0d0) then
            fvec(i) = ktc(i)
         else
            fvec(i) = 0d0
         endif
      end do
      alpha_out = alpha
      return
      contains
      ! construct G as sum_k k0(:,k) k0(:,k)^T for columns where (K0^T*c)(k)>0
      subroutine compute_gmat(k0,m,n,c,g)
         integer,intent(in) :: m,n
         double precision,intent(in) :: k0(m,n),c(m)
         double precision,intent(out) :: g(m,m)
         double precision,dimension(n) :: prod
         double precision,dimension(m,n) :: kpos
         integer :: l,poscols
         prod = matmul(transpose(k0),c)
         poscols = 0
         do l=1,n
            if (prod(l) > 0d0) then
               poscols = poscols + 1
               kpos(:,poscols) = k0(:,l)
            end if
         end do
         if (poscols > 0) then
            ! sum of outer products -> matmul of selected columns
            g = matmul(kpos(:,1:poscols), transpose(kpos(:,1:poscols)))
         else
            g = 0d0
         end if
      end subroutine compute_gmat
      subroutine compute_grad(g,c,mvec,alpha,grad)
         double precision,intent(in) :: g(m,m),c(m),mvec(m),alpha
         double precision,intent(out) :: grad(m)
         grad = matmul(g,c) + alpha*c - mvec
      end subroutine compute_grad
      subroutine compute_hessian(g,alpha,h)
         double precision,intent(in) :: g(m,m),alpha
         double precision,intent(out) :: h(m,m)
         h = g
         do i=1,m
            h(i,i) = h(i,i) + alpha
         end do
      end subroutine compute_hessian
      end subroutine venk_brd
