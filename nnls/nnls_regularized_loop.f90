      subroutine nnls_regularized_loop(a,m,n,q,b,x,rnorm,w,zz,idx, &
              mode,maxiter,lambda) 
          !f2py threadsafe
          integer, intent(in):: m, n, q, maxiter
          ! we no longer need to specify a or b as intent copy,
          ! since we are going to manually copy them
          double precision, intent(in):: a(m,n),b(m,q),lambda
          double precision, intent(out):: x(n,q)
          double precision x_temp(n)
          double precision w(*), zz(*)
          double precision a_prime(m+n,n), b_prime(m+n)
          double precision b_prime_(m+n,q)
          double precision, intent(out):: rnorm(m)
          double precision rnorm_temp
          integer idx(*), j, k, p
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
          b_prime_(1:m,1:q) = b(1:m,1:q)
          b_prime_(m+1:m+n,1:q) = 0d0
          do p=1,q
              b_prime(1:m+n) = b_prime_(1:m+n,p)
              call nnls(a_prime,m+n,n,b_prime,x_temp,rnorm_temp,w,zz,idx,mode,maxiter)
              x(1:n,p) = x_temp(1:n)
              rnorm(p) = rnorm_temp
          end do
      end subroutine
