  subroutine venk_brd(initial_alpha,k0_mat,m,n,m_r,f,alpha_out,tol,maxiter)
    !f2py threadsafe
    ! implement the Butler-Reeds-Dawson solver in Fortran
    integer,intent(in) :: m,n,maxiter
    double precision,intent(in) :: initial_alpha,tol
    double precision,intent(in) :: k0_mat(m,n),m_r(m)
    double precision,intent(out) :: f(n),alpha_out
    double precision,allocatable :: c(:)
    ! goes to sub
    double precision,allocatable :: tempvec(:)
    double precision :: alpha,alpha_new,sqrt_n
    integer :: iter,j
    double precision :: norm_mr
    external dgemv
    sqrt_n = sqrt(dble(m))
    alpha = initial_alpha

    allocate(c(m),tempvec(n))
    c = 1.0d0
    norm_mr = norm2(m_r)
    do iter=1,maxiter
      call venk_nnls(k0_mat,m_r,c,alpha,m,n)
      alpha_new = sqrt_n / norm2(c)
      write(*,*) 'alpha iteration', iter, 'value', alpha_new
      if (abs(alpha_new-alpha)/alpha < tol) exit
      alpha = alpha_new
    end do
    call dgemv('T',m,n,1.0d0,k0_mat,m,c,1,0.0d0,tempvec,1)
    do j=1,n
      if (tempvec(j) > 0.d0) then
        f(j) = tempvec(j)
      else
        f(j) = 0.d0
      end if
    end do
    deallocate(c,tempvec)
    alpha_out = alpha
    return
  end subroutine venk_brd
  subroutine venk_nnls(k0_mat,m_r,c,alpha,m,n)
    integer,intent(in)::m,n
    integer :: info, k
    double precision,allocatable :: tempvec(:)
    double precision,intent(in) :: alpha,k0_mat(m,n),m_r(m)
    double precision,intent(inout) :: c(m)
    double precision,allocatable :: g_mat(:,:),c_new(:),grad(:),newgrad(:)
    double precision,allocatable :: h_mat(:,:),h_mat_copy(:,:),k0_t(:)
    integer,allocatable :: piv(:)
    double precision,allocatable :: delta_c(:)
    double precision :: chi_old,chi_new,s,denom,norm_grad,norm_mr
    external dgesv, dgemm
    allocate(delta_c(m), tempvec(m))
    allocate(g_mat(m,m),c_new(m),grad(m),newgrad(m),h_mat(m,m),h_mat_copy(m,m),k0_t(n),piv(m))
    norm_mr = norm2(m_r)
    !write(*,*) 'venk_nnls called with', alpha, 'and initial c'
    !do j=1,m
    !  write(*,*) c(j)
    !end do
    do j=1,500
      ! IN: k0_mat,m,n,c OUT: g_mat,k0_t
      call compute_g(k0_mat,m,n,c,g_mat,k0_t)
      ! IN: g_mat,c,m,alpha,m_r OUT: grad
      call compute_grad(g_mat,c,m,alpha,m_r,grad)
      ! IN: g_mat,m,alpha OUT: h_mat
      call add_diag(g_mat,m,alpha,h_mat)
      delta_c = grad
      h_mat_copy = h_mat
      ! IN: m, 1, m, m INOUT: h_mat,delta_c OUT: piv,info
      call dgesv(m,1,h_mat_copy,m,piv,delta_c,m,info)
      call dgemm('N','N',m,1,m,1.0d0,h_mat,m,delta_c,m,0.0d0,tempvec,m)
      denom = dot_product(delta_c,tempvec)
      if (denom == 0.d0) denom = 1.d-12
      s = dot_product(delta_c,grad) / denom
      c_new = c - s*delta_c
      ! IN: g_mat,c_new,m,alpha,m_r OUT: newgrad
      call compute_grad(g_mat,c_new,m,alpha,m_r,newgrad)
      ! IN: c,g_mat,alpha,m_r OUT: chi_old
      call chi_func(c,g_mat,alpha,m_r,chi_old)
      ! IN: c_new,g_mat,alpha,m_r OUT: chi_new
      call chi_func(c_new,g_mat,alpha,m_r,chi_new)
      k=0
      do while (chi_new >= chi_old .and. k < 20)
        k = k + 1
        c_new = c - s*(0.5d0**k)*delta_c
        ! IN: g_mat,c_new,m,alpha,m_r OUT: newgrad
        call compute_grad(g_mat,c_new,m,alpha,m_r,newgrad)
        ! IN: c_new,g_mat,alpha,m_r OUT: chi_new
        call chi_func(c_new,g_mat,alpha,m_r,chi_new)
      end do
      c = c_new
      grad = newgrad
      norm_grad = norm2(grad)
      if (norm_grad/norm_mr < 1.d-8) exit
    end do
    deallocate(g_mat,c_new,newgrad,h_mat,h_mat_copy,k0_t,piv,delta_c,tempvec)
    contains
      subroutine chi_func(c_forchi,g_mat_forchi,alpha_forchi,m_r_forchi,val_forchi)
        double precision,intent(in)::c_forchi(:),g_mat_forchi(size(c_forchi),size(c_forchi)),alpha_forchi,m_r_forchi(:)
        double precision,intent(out)::val_forchi
        double precision :: temp(size(c_forchi))
        call dgemm('N','N',size(c_forchi),1,size(c_forchi),1.0d0, &
            g_mat_forchi,size(c_forchi),c_forchi,size(c_forchi),0.0d0, &
            temp,size(c_forchi))
        val_forchi = dot_product(c_forchi,temp) &
          & + alpha_forchi*dot_product(c_forchi,c_forchi) &
          & - dot_product(c_forchi,m_r_forchi)
      end subroutine chi_func
  end subroutine venk_nnls
  subroutine compute_g(k0_mat,m,n,c,g_mat,k0_t)
    integer,intent(in)::m,n
    double precision,intent(in)::k0_mat(m,n),c(m)
    double precision,intent(out)::g_mat(m,m),k0_t(n)
    integer :: p
    g_mat = 0.d0
    k0_t = 0.d0
    do p=1,n
      k0_t(p) = dot_product(k0_mat(:,p),c)
      if (k0_t(p) > 0.d0) then
        ! IN: k0_mat(:,p) INOUT: g_mat
        call add_outer(k0_mat(:,p),g_mat)
      end if
    end do
    contains
      subroutine add_outer(col,a_mat)
        double precision,intent(in)::col(:)
        double precision,intent(inout)::a_mat(:,:)
        integer :: i,j
        do i=1,size(col)
          do j=1,size(col)
            a_mat(i,j) = a_mat(i,j) + col(i)*col(j)
          end do
        end do
      end subroutine add_outer
  end subroutine compute_g
  subroutine compute_grad(g_mat,c,m,alpha,m_r,grad)
    integer,intent(in)::m
    double precision,intent(in)::g_mat(m,m),c(m),alpha,m_r(m)
    double precision,intent(out)::grad(m)
    double precision :: temp(m)
    call dgemm('N','N',m,1,m,1.0d0,g_mat,m,c,m,0.0d0,temp,m)
    grad = temp + alpha*c - m_r
  end subroutine compute_grad
  subroutine add_diag(g_mat,m,alpha,h_mat)
    integer,intent(in)::m
    double precision,intent(in)::g_mat(m,m),alpha
    double precision,intent(out)::h_mat(m,m)
    integer::i
    h_mat = g_mat
    do i=1,m
      h_mat(i,i) = h_mat(i,i) + alpha
    end do
  end subroutine add_diag
