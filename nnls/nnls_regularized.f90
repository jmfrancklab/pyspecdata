      subroutine nnls_regularized(a,m,n,b,x,rnorm,w,zz,idx, &
              mode,maxiter,lambda) 
          !f2py threadsafe
          integer, intent(in):: m, n, maxiter
          ! we no longer need to specify a or b as intent copy,
          ! since we are going to manually copy them
          double precision, intent(in):: a(m,n),b(*),lambda
          double precision, intent(out):: x(n)
          double precision w(*), zz(*)
          double precision a_prime(m+n,n), b_prime(m+n)
          double precision, intent(out):: rnorm
          integer idx(*), j, k
          integer, intent(out):: mode
          a_prime(1:m,1:n) = a
          do j=1,n
              do k=1,n
                  if (j == k) then
                      a_prime(m+j,k) = lambda
                  else
                      a_prime(m+j,k) = 0.0d0
                  endif
              end do
          end do
          b_prime(1:m) = b(1:m)
          b_prime(m+1:m+n) = 0d0
          call nnls(a_prime,m+n,n,b_prime,x,rnorm,w,zz,idx,mode,maxiter)
      end subroutine
