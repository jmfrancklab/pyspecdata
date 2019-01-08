      subroutine nnls_regularized_loop(a,m,n,b,x,rnorm,w,zz,idx, &
              mode,maxiter,lambda,l) 
          integer, intent(in):: m, n, maxiter, l
          ! we no longer need to specify a or b as intent copy,
          ! since we are going to manually copy them
          double precision, intent(in):: a(m,n),b(*),lambda(l)
          double precision, intent(out):: x(l,n)
          double precision x_temp(n)
          double precision w(*), zz(*)
          double precision a_prime(m+n,n), b_prime(m+n)
          double precision, intent(out):: rnorm(l)
          double precision rnorm_temp
          integer idx(*), i, j, k
          integer, intent(out):: mode
          a_prime(1:m,1:n) = a
          do i=0,l
              do j=1,n
                  do k=1,n
                      if (j == k) then
                          a_prime(m+j,k) = lambda(i)
                      else
                          a_prime(m+j,k) = 0.0d0
                      endif
                  end do
              end do
              b_prime(1:m) = b(1:m)
              b_prime(m+1:m+n) = 0d0
              call nnls(a_prime,m+n,n,b_prime,x_temp,rnorm_temp,w,zz,idx,mode,maxiter)
              x(i,1:n) = x_temp(1:n)
              rnorm(i) = rnorm_temp
          end do
      end subroutine
