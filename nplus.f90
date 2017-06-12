!****************************************************************************
! module: prec
! COMMENT: Define parameter for NRG subroutine
!****************************************************************************

    module prec
    
        implicit none
        integer,parameter::Rkind=8
        
    end module prec

!****************************************************************************
! module: precnpl
! COMMENT: Define parameter for NRG subroutine
!****************************************************************************

    module precnpl

        use prec
        implicit none

        integer,parameter:: Nit=200, Nnit=30         ! Max iteration numbers.
        integer,parameter:: Nb0=2000, Ns=80, Nb=30, Dim=2400

        real(kind=Rkind), parameter:: Pi=3.141526535_Rkind
        real(kind=Rkind), parameter:: Lamda=10.0_Rkind
        real(kind=Rkind), parameter:: wc=1.0_Rkind, s=0.5_Rkind
        
        real(kind=Rkind):: dta_loc=0.1_Rkind, epi_loc=0.0_Rkind   !  local term
        real(kind=Rkind):: g1_cps=0.0_Rkind, g2_cps=1.0_Rkind      !  coupling strenth         
        
    end module precnpl

!****************************************************************************
! PROGRAM: NPL_NRG
! TYPE   : Main
! PURPOSE: Nonlear + Linear 
!****************************************************************************
    program NPL_NRG
    use precnpl
    implicit none

    integer:: i, j, k, ierror, ierror2, jishu

    real(kind=Rkind):: time_begin, time_end
    real(kind=Rkind):: eta0, epi0, t0, epin, tn
    real(kind=Rkind):: alpha
    real(kind=Rkind),dimension(0:Nit):: epi, t

    real(kind=Rkind),dimension(Ns):: Eigenv
    real(kind=Rkind),dimension(Ns, Ns):: bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix
    real(kind=Rkind),dimension(50, 0:Nit, Ns):: Eigenv_all, w_all, cw_all, cwx_all, cwsx_all

!--- Begin ---
    call cpu_time(time_begin)
    alpha=0.0000265392_Rkind

!--- Begin Caculation	    
    Eigenv_all=0.0_Rkind
    w_all=0.0_Rkind
    cw_all=0.0_Rkind
    cwx_all=0.0_Rkind
    cwsx_all=0.0_Rkind

!$omp do
    do j=1, 1

        call param(alpha, eta0, epi, t)
        
        !--- First iteration
        epi0=epi(0)
        call diag_H0(eta0, epi0, Eigenv, bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix)
        print *, "<qz>:", qzmatrix(1, 1), qxmatrix(1, 1), abs(b0matrixa(1, 1)+b0matrixc(1, 1)), x2matrix(1, 1), nmatrix(1, 1)

        do k=1, Ns
            Eigenv_all(j, 0, k)=Eigenv(k)                  !--- Save eigenv
            w_all(j, 0, k)=Eigenv(k)
            cw_all(j, 0, k)=0.5_Rkind*(qzmatrix(1, k)**2)
            cwsx_all(j, 0, k)=0.5_Rkind*(qxmatrix(1, k)**2)
            cwx_all(j, 0, k)=0.5_Rkind*((b0matrixa(1, k)+b0matrixc(1, k))**2)
        end do

        !--- Later iterations
        do i=1, Nnit
        
            epin=epi(i)
            tn=Lamda*t(i-1)
            call diag_Hn(epin, tn, Eigenv, bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix)
            print *, "<qz>:", qzmatrix(1, 1), abs(b0matrixa(1, 1)+b0matrixc(1, 1)), x2matrix(1, 1), nmatrix(1, 1)

            do k=1, Ns
                Eigenv_all(j, i, k)=Eigenv(k)                   !--- Save eigenv
                w_all(j, i, k)=Eigenv(k)*(Lamda**(-i))
                cw_all(j, i, k)=0.5_Rkind*(qzmatrix(1, k)**2)
                cwsx_all(j, i, k)=0.5_Rkind*(qxmatrix(1, k)**2)
                cwx_all(j, i, k)=0.5_Rkind*((b0matrixa(1, k)+b0matrixc(1, k))**2)
            end do

        end do
        
        !--- Physical information. Saved to disk.
        OPEN (UNIT=2, FILE='m.TXT', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror2)
        WRITE (2,900) alpha, qzmatrix(1, 1), qxmatrix(1, 1), (b0matrixa(1, 1)+b0matrixc(1, 1)), x2matrix(1, 1), nmatrix(1, 1), (abs(b0matrixa(1, 1)+b0matrixc(1, 1))**2)/x2matrix(1, 1)
        900 FORMAT (F, F, F, F, F, F, F)

        alpha = alpha + 0.0000001_Rkind
    
    end do
!$omp end do

!--- Flow diagrams. Saved to disk.    
    jishu=0
    do i=0, Nnit

        OPEN (UNIT=1, FILE="flow.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
        WRITE (1,700) jishu, Eigenv_all(1, i, 2), Eigenv_all(2, i, 2), Eigenv_all(3, i, 2), Eigenv_all(4, i, 2), Eigenv_all(5, i, 2), Eigenv_all(6, i, 2), Eigenv_all(7, i, 2), Eigenv_all(8, i, 2)
        700 FORMAT (I, F, F, F, F, F, F, F, F)

        jishu=jishu+1
    
    end do

    !call calc_cw(w_all, cw_all)
    !call calc_cwx(w_all, cwx_all)
    call calc_cwsx(w_all, cwsx_all)

!--- End----
    call cpu_time(time_end)
    print *, 'Time of operation was : ', time_end - time_begin, ' seconds' 

    end program NPL_NRG


!========+=========+===============+=========+=$
! PROGRAM: calc_cwsx
! TYPE   : subroutine
! PURPOSE: 
!========+=========+===============+=========+=$    
    subroutine calc_cwsx(w_all, cwsx_all)
    use precnpl
    implicit none
    
    integer:: i, j, k, ierror3, ierror4, ierror8
    real(kind=Rkind):: sum, w, res_cwsx, dlt_tmp
    
    real(kind=Rkind),dimension(50, 0:Nit, Ns):: w_all, cwsx_all
        
    do i=0, Nnit
        print *, 'w_max: ', w_all(1, i, Ns)
    end do

    do i=0, Nnit
        do j=1, Ns

        OPEN (UNIT=18, FILE="cwsx_origin.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror8)
        WRITE (18,2018)  w_all(1, i, j), cwsx_all(1, i, j)
        2018 FORMAT (F, F)
        
        end do
    end do

    do i=1, Nnit
        do j=0, i-1
        
            k=1
            do
            if(w_all(1, j, k)>w_all(1, i, Ns)) exit
            cwsx_all(1, j, k)=cwsx_all(1, j, k)*(w_all(1, j, k)/w_all(1, i, Ns))
            k=k+1
            end do
            
            
        end do
        
        do j=1, Ns
            cwsx_all(1, i, j)=cwsx_all(1, i, j)*(1-(w_all(1, i, j)/w_all(1, i, Ns)))
        end do
    end do

    do i=0, Nnit
        print *, 'cw_0: ', cwsx_all(1, i, 1)
    end do
    
    sum=0.0_Rkind
    do i=0, Nnit
        do j=1, Ns
        sum=sum+cwsx_all(1, i, j)
        end do
    end do
    print *, 'sum: ', sum
    
    w=0.0_Rkind
    do k=1, 90
    
        res_cwsx=0.0_Rkind
        do i=0, Nnit
            do j=1, Ns
            
                if(w_all(1, i, j)==0.0_Rkind) then
                
                    dlt_tmp=0.0_Rkind
                    
                else    
                
                    dlt_tmp=0.6_Rkind**2
                    dlt_tmp=exp(-(dlt_tmp/4))/(0.6_Rkind*w_all(1, i, j)*sqrt(Pi))
                    dlt_tmp=dlt_tmp*exp(-((log(w)-log(w_all(1, i, j)))**2)/(0.6_Rkind**2))      
                    
                end if  
                
                res_cwsx=res_cwsx+cwsx_all(1, i, j)*dlt_tmp
            end do
        end do
    
    
            OPEN (UNIT=28, FILE="res_cwsx.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror3)
            WRITE (28,2028)  w, res_cwsx
            2028 FORMAT (F, F)
    
        w=0.00000000000001_Rkind*1.5_Rkind**k
    end do    
    
    do i=0, Nnit
        do j=1, Ns

        OPEN (UNIT=29, FILE="cwsx.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror4)
        WRITE (29 ,2029)  w_all(1, i, j), cwsx_all(1, i, j)
        2029 FORMAT (F, F)
        
        if(cwsx_all(1, i, j)>0.01_Rkind) then
        print *, 'i  j', i, j, w_all(1, i, j), cwsx_all(1, i, j)
        end if

        end do
    end do

    end subroutine calc_cwsx

!========+=========+===============+=========+=$
! PROGRAM: calc_cw
! TYPE   : subroutine
! PURPOSE: 
!========+=========+===============+=========+=$    
    subroutine calc_cw(w_all, cw_all)
    use precnpl
    implicit none
    
    integer:: i, j, k, ierror3, ierror4, ierror8
    real(kind=Rkind):: sum, w, res_cw, dlt_tmp
    
    real(kind=Rkind),dimension(50, 0:Nit, Ns):: w_all, cw_all
        
    do i=0, Nnit
        print *, 'w_max: ', w_all(1, i, Ns)
    end do

    do i=0, Nnit
        do j=1, Ns

        OPEN (UNIT=8, FILE="cw_origin.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror8)
        WRITE (8,2000)  w_all(1, i, j), cw_all(1, i, j)
        2000 FORMAT (F, F)
        
        end do
    end do

    do i=1, Nnit
        do j=0, i-1
        
            k=1
            do
            if(w_all(1, j, k)>w_all(1, i, Ns)) exit
            cw_all(1, j, k)=cw_all(1, j, k)*(w_all(1, j, k)/w_all(1, i, Ns))
            k=k+1
            end do
            
            
        end do
        
        do j=1, Ns
            cw_all(1, i, j)=cw_all(1, i, j)*(1-(w_all(1, i, j)/w_all(1, i, Ns)))
        end do
    end do

    do i=0, Nnit
        print *, 'cw_0: ', cw_all(1, i, 1)
    end do
    
    sum=0.0_Rkind
    do i=0, Nnit
        do j=1, Ns
        sum=sum+cw_all(1, i, j)
        end do
    end do
    print *, 'sum: ', sum
    
    w=0.0_Rkind
    do k=1, 90
    
        res_cw=0.0_Rkind
        do i=0, Nnit
            do j=1, Ns
            
                if(w_all(1, i, j)==0.0_Rkind) then
                
                    dlt_tmp=0.0_Rkind
                    
                else    
                
                    dlt_tmp=0.6_Rkind**2
                    dlt_tmp=exp(-(dlt_tmp/4))/(0.6_Rkind*w_all(1, i, j)*sqrt(Pi))
                    dlt_tmp=dlt_tmp*exp(-((log(w)-log(w_all(1, i, j)))**2)/(0.6_Rkind**2))      
                    
                end if  
                
                res_cw=res_cw+cw_all(1, i, j)*dlt_tmp
            end do
        end do
    
    
            OPEN (UNIT=3, FILE="res_cw.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror3)
            WRITE (3,900)  w, res_cw
            900 FORMAT (F, F)
    
        w=0.00000000000001_Rkind*1.5_Rkind**k
    end do    
    
    do i=0, Nnit
        do j=1, Ns

        OPEN (UNIT=4, FILE="cw.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror4)
        WRITE (4,1000)  w_all(1, i, j), cw_all(1, i, j)
        1000 FORMAT (F, F)
        
        if(cw_all(1, i, j)>0.01_Rkind) then
        print *, 'i  j', i, j, w_all(1, i, j), cw_all(1, i, j)
        end if

        end do
    end do

    end subroutine calc_cw
    
    
!========+=========+===============+=========+=$
! PROGRAM: calc_cwx
! TYPE   : subroutine
! PURPOSE: 
!========+=========+===============+=========+=$    
    subroutine calc_cwx(w_all, cwx_all)
    use precnpl
    implicit none
    
    integer:: i, j, k, ierror5, ierror6
    real(kind=Rkind):: sum, w, res_cwx, dlt_tmp
    
    real(kind=Rkind),dimension(50, 0:Nit, Ns):: w_all, cwx_all
    
    
    do i=0, Nnit
        print *, 'w_max: ', w_all(1, i, Ns)
    end do
    
    do i=1, Nnit
        do j=0, i-1
        
            k=1
            do
            if(w_all(1, j, k)>w_all(1, i, Ns)) exit
            cwx_all(1, j, k)=cwx_all(1, j, k)*(w_all(1, j, k)/w_all(1, i, Ns))
            k=k+1
            end do
            
            
        end do
        
        do j=1, Ns
            cwx_all(1, i, j)=cwx_all(1, i, j)*(1-(w_all(1, i, j)/w_all(1, i, Ns)))
        end do
    end do

    do i=0, Nnit
        print *, 'cwx_0: ', cwx_all(1, i, 1)
    end do
    
    sum=0.0_Rkind
    do i=0, Nnit
        do j=1, Ns
        sum=sum+cwx_all(1, i, j)
        end do
    end do
    print *, 'sum: ', sum
    
    do k=1, 90
    w=0.00000000000001_Rkind*1.5_Rkind**k
    
        res_cwx=0.0_Rkind
        do i=0, Nnit
            do j=1, Ns
            
                if(w_all(1, i, j)<0.0000000000000001_Rkind) then
                
!                   dlt_tmp=(0.5_Rkind/Pi)*0.6_Rkind/((w-w_all(1, i, j))**2+0.6_Rkind**2)
                    dlt_tmp=0.0_Rkind
                    
                else    
                
                    dlt_tmp=0.6_Rkind**2
                    dlt_tmp=exp(-(dlt_tmp/4))/(0.6_Rkind*w_all(1, i, j)*sqrt(Pi))
                    dlt_tmp=dlt_tmp*exp(-((log(w)-log(w_all(1, i, j)))**2)/(0.6_Rkind**2))      
                    
                end if  
                
                res_cwx=res_cwx+cwx_all(1, i, j)*dlt_tmp
            end do
        end do
    
    
            OPEN (UNIT=5, FILE="res_cwx.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror5)
            WRITE (5,1100)  w, res_cwx
            1100 FORMAT (F, F)
    
    end do
    
    
    
    
    do i=0, Nnit
        do j=1, Ns

        OPEN (UNIT=6, FILE="cwx.TXT", STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror6)
        WRITE (6,1200)  w_all(1, i, j), cwx_all(1, i, j)
        1200 FORMAT (F, F)
        
        if(cwx_all(1, i, j)>0.01_Rkind) then
        print *, 'i  j', i, j, w_all(1, i, j), cwx_all(1, i, j)
        end if

        end do
    end do

    end subroutine calc_cwx

!========+=========+===============+=========+=$
! PROGRAM: param
! TYPE   : subroutine
! PURPOSE: Produce NRG Paranmeters
!========+=========+===============+=========+=$	
    subroutine param(alpha, eta0, epi, t)
    use precnpl
    implicit none
    
!--- Define parameter----
    integer:: i, j, k, Isign

    real(kind=Rkind):: alpha
    real(kind=Rkind):: eta0, tmp1, tmp2
    real(kind=Rkind),dimension(0:Nit):: xi, gama2, epi, t, orth
    real(kind=Rkind),dimension(0:Nit,0:Nit):: Umatrix

    do i=0, Nit-1
        gama2(i)=(2*Pi*alpha)*(wc**2)*(1-lamda**(-(s+1)))*(lamda**(-i*(s+1)))/(s+1)
        xi(i)=((s+1)/(s+2))*((1-lamda**(-(s+2)))/(1-lamda**(-(s+1))))*wc*(lamda**(-i))
    end do

    eta0=0.0_Rkind
    do i=0, Nit-1
        eta0=eta0+gama2(i)
    end do


    Umatrix=0.0_Rkind
    epi=0.0_Rkind
    do i=0, Nit-1
        tmp1=gama2(i)/eta0
        Umatrix(0, i)=sqrt(tmp1)
        epi(0)=epi(0)+xi(i)*tmp1
    end do

    t=0.0_Rkind
    do i=0, Nit-1
        t(0)=((xi(i)-epi(0))*Umatrix(0, i))**2+t(0)
    end do
    t(0)=sqrt(t(0))

    do i=0, Nit-1
        Umatrix(1, i)=(xi(i)-epi(0))*Umatrix(0, i)/t(0)
        epi(1)=epi(1)+xi(i)*(Umatrix(1, i)**2)
    end do

    do i=0, Nit-1
        t(1)=((xi(i)-epi(1))*Umatrix(1, i)-t(0)*Umatrix(0, i))**2+t(1)
    end do
    t(1)=sqrt(t(1))

!----------iteration-----------
    do i=2, Nit-1

        do j=0, Nit-1
            Umatrix(i, j)=((xi(j)-epi(i-1))*Umatrix(i-1, j)-t(i-2)*Umatrix(i-2, j))/t(i-1)
        end do

!----------orth-----------
        orth=0.0_Rkind
        do j=0, i-1

            tmp2=0.0_Rkind
            do k=0, Nit-1
                tmp2=tmp2+Umatrix(i, k)*Umatrix(j, k)
            end do

            do k=0, Nit-1
                orth(k)=orth(k)+tmp2*Umatrix(j, k)
            end do

        end do

        do j=0, Nit-1
            Umatrix(i, j)=Umatrix(i, j)-orth(j)
        end do

        tmp2=0.0_Rkind
        do j=0, Nit-1
            tmp2=tmp2+Umatrix(i, j)*Umatrix(i, j)
        end do
        
        tmp2=sqrt(tmp2)
        do j=0, Nit-1
            Umatrix(i, j)=Umatrix(i, j)/tmp2
        end do

        do j=0, Nit-1
            epi(i)=epi(i)+xi(j)*(Umatrix(i, j)**2)
        end do

        do j=0, Nit-1
            t(i)=((xi(j)-epi(i))*Umatrix(i, j)-t(i-1)*Umatrix(i-1, j))**2+t(i)
        end do
        t(i)=sqrt(t(i))

    end do

    do i=0, Nit-1
        t(i)=t(i)*(lamda**i)
        epi(i)=epi(i)*(lamda**i)
        print *, i, t(i), epi(i)
    end do

    Isign=15
!        print *, Isign-4, t(Isign-4)
!        print *, Isign-3, t(Isign-3)
!        print *, Isign-2, t(Isign-2)
!        print *, Isign-1, t(Isign-1)
!        print *, Isign, t(Isign)
!        print *, Isign+1, t(Isign+1)
!        print *, Isign+2, t(Isign+2)
!        pause

    do i=Isign+1, Nit-1
        epi(i)=epi(Isign)
        t(i)=t(Isign)
    end do

    end subroutine param

!========+=========+===============+=========+=$
! PROGRAM: diag_H0
! TYPE   : subroutine
! PURPOSE: 
!========+=========+===============+=========+=$	
    subroutine diag_H0(eta0, epi0, Eigenv, bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix)
    use precnpl
    implicit none

    integer:: i, j, k

    real(kind=Rkind):: eta0, epi0, t0, Coeff0, Coeff1, tmp, sum1, sum2

    real(kind=Rkind),dimension(Ns):: Eigenv
    real(kind=Rkind),dimension(Nb0):: vec_mid, vec_mid1, vec_mid2
    real(kind=Rkind),dimension(Ns, Ns):: bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix
    real(kind=Rkind),dimension(Nb0, Nb0):: H0
    real(kind=Rkind),dimension(Nb0/2, Nb0/2):: H01

    integer:: cevec=1
    integer:: ierr
    real(kind=Rkind),dimension(Nb0):: eigvalue, fv1, fv2
    real(kind=Rkind),dimension(Nb0, Nb0):: eigvector

!    print *, "eta0, epi0:", eta0, epi0
    H0=0.0_Rkind

!----------loc term-----------
    do i=1, Nb0/2
        H0(i, i)=H0(i, i)+(epi_loc/2)
        H0(i+Nb0/2, i+Nb0/2)=H0(i+Nb0/2, i+Nb0/2)-(epi_loc/2)
        H0(i+Nb0/2, i)=H0(i+Nb0/2, i)-(dta_loc/2)
        H0(i, i+Nb0/2)=H0(i, i+Nb0/2)-(dta_loc/2)
    end do

!----------linear term-----------
    Coeff0=g1_cps*sqrt(eta0/Pi)/2

    do i=1, Nb0/2-1
        tmp=i
        tmp=sqrt(tmp)
        tmp=Coeff0*tmp
        H0(i, i+1)=H0(i, i+1)+tmp
        H0(i+1, i)=H0(i+1, i)+tmp
        H0(i+Nb0/2, i+Nb0/2+1)=H0(i+Nb0/2, i+Nb0/2+1)-tmp
        H0(i+Nb0/2+1, i+Nb0/2)=H0(i+Nb0/2+1, i+Nb0/2)-tmp
    end do    
        
    
!----------nonlinearg term-----------
    Coeff1=g2_cps*eta0/Pi/2

    H01=0.0_Rkind
    do i=1, Nb0/2-1
        tmp=i
        tmp=sqrt(tmp)
        H01(i, i+1)=H01(i, i+1)+tmp
        H01(i+1, i)=H01(i+1, i)+tmp
    end do

    do i=1, Nb0/2
        do j=1, Nb0/2
            tmp=0.0_Rkind
            do k=1, Nb0/2
                tmp=tmp+H01(i, k)*H01(k, j)
            end do
            H0(i, j)=H0(i, j)+Coeff1*tmp
            H0(i+Nb0/2, j+Nb0/2)=H0(i+Nb0/2, j+Nb0/2)-Coeff1*tmp
        end do
    end do

!----------bath term-----------
    do i=1, Nb0/2
        H0(i, i)=H0(i, i)+epi0*(i-1)
        H0(i+Nb0/2, i+Nb0/2)=H0(i+Nb0/2, i+Nb0/2)+epi0*(i-1)
    end do

    do i=1, Nb0
        do j=1, Nb0
        if(H0(i, j)/=H0(j, i)) then
        print *, "e"
        pause
        end if
        end do
    end do

!    print *, H0

!-------------- diagonalize Hamilt matrix ------------------
!---- Call RS subroutine --------
!--using slatec in 1976: slower, but independent of compiler ---- 
!          call  rs (Nb0, Nb0, H0, eigvalue, &
!                            cevec, eigvector, fv1, fv2, ierr)
!--using Lapack: faster, but must with intel + Lapack ----
          call  rs_cpu (Nb0, Nb0, H0, eigvalue, &
                    cevec, eigvector, fv1, fv2, ierr)
!---------------------
    if (ierr/=0) then
        print *, 'error'
        stop
    endif
!----------------------------------------------------

!---- Output eigenvalues --------    
    do i=1, Ns
        Eigenv(i)=eigvalue(i)
    end do

    tmp=Eigenv(1)
    do i=1, Ns
        Eigenv(i)=Eigenv(i)-tmp
    end do

!---- Write b0 matrix --------     
    do i=1, Ns

        vec_mid1=0.0_Rkind
        do j=1, Nb0/2-1
            tmp=j
            tmp=sqrt(tmp)
            vec_mid1(j)=tmp*eigvector(j+1, i)
            vec_mid1(j+Nb0/2)=tmp*eigvector(j+Nb0/2+1, i)
        end do
        
        do j=1, Ns

            sum1=0.0_Rkind
            do k=1, Nb0
                sum1=sum1+vec_mid1(k)*eigvector(k, j)
            end do

            bmatrixa(j, i)=sum1
        end do
    end do

!---- Write b0+ matrix --------    
    do i=1, Ns
        do j=1, Ns
            bmatrixc(i, j)=bmatrixa(j, i)
            b0matrixa(j, i)=bmatrixa(j, i)
            b0matrixc(i, j)=bmatrixa(j, i)
        end do
    end do

!---- Write x^2 matrix --------     
    do i=1, Ns

        vec_mid=0.0_Rkind
        vec_mid1=0.0_Rkind
        vec_mid2=0.0_Rkind
        do j=1, Nb0/2-1
            tmp=j
            tmp=sqrt(tmp)
            vec_mid1(j)=tmp*eigvector(j+1, i)
            vec_mid1(j+Nb0/2)=tmp*eigvector(j+Nb0/2+1, i)
            vec_mid2(j+1)=tmp*eigvector(j, i)
            vec_mid2(j+Nb0/2+1)=tmp*eigvector(j+Nb0/2, i)
        end do

        do j=1, Nb0
            vec_mid(j)=vec_mid1(j)+vec_mid2(j)
        end do        

        vec_mid1=0.0_Rkind
        vec_mid2=0.0_Rkind
        do j=1, Nb0/2-1
            tmp=j
            tmp=sqrt(tmp)
            vec_mid1(j)=tmp*vec_mid(j+1)
            vec_mid1(j+Nb0/2)=tmp*vec_mid(j+Nb0/2+1)
            vec_mid2(j+1)=tmp*vec_mid(j)
            vec_mid2(j+Nb0/2+1)=tmp*vec_mid(j+Nb0/2)
        end do

        do j=1, Nb0
            vec_mid(j)=vec_mid1(j)+vec_mid2(j)
        end do

        do j=1, Ns

            sum1=0.0_Rkind
            do k=1, Nb0
                sum1=sum1+vec_mid(k)*eigvector(k, j)
            end do

            x2matrix(j, i)=sum1
        end do
    end do


!---- Write <a+a> matrix --------     
    do i=1, Ns

        do j=1, Ns

            sum1=0.0_Rkind
            do k=1, Nb0/2
                tmp=k-1
                sum1=sum1+tmp*eigvector(k, i)*eigvector(k, j)+tmp*eigvector(k+Nb0/2, i)*eigvector(k+Nb0/2, j)
            end do

            nmatrix(j, i)=sum1
        end do
    end do

!---- Write qx matrix --------    
    do i=1, Ns
        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Nb0/2
                sum2=sum2+eigvector(k+Nb0/2, i)*eigvector(k, j)+eigvector(Nb0/2, i)*eigvector(k+Nb0/2, j)
            end do
            qxmatrix(j, i)=sum2

        end do
    end do

!---- Write qz matrix --------    
    do i=1, Ns
        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Nb0/2
                sum2=sum2+eigvector(k, i)*eigvector(k, j)+eigvector(k+Nb0/2, i)*(-eigvector(k+Nb0/2, j))
            end do
            qzmatrix(j, i)=sum2

        end do
    end do

!---- Test Hermite --------    
!    do i=1, Ns
!        do j=1, Ns
!            if(qzmatrix(i, j)/=qzmatrix(j, i)) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

!    do i=1, Ns
!        do j=1, Ns
!            if(b0matrixa(i, j)/=b0matrixc(j, i)) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

    end subroutine diag_H0

!========+=========+===============+=========+=$
! PROGRAM: diag_Hn
! TYPE   : subroutine
! PURPOSE: 
!========+=========+===============+=========+=$	
    subroutine diag_Hn(epin, tn, Eigenv, bmatrixc, bmatrixa, qzmatrix, qxmatrix, b0matrixa, b0matrixc, x2matrix, nmatrix)
    use precnpl
    implicit none

    integer:: i, j, k, l
    real(kind=Rkind):: epin, tn, tmp, sum1, sum2

    real(kind=Rkind),dimension(Nb, Nb):: bmatrixan1, bmatrixcn1

    real(kind=Rkind),dimension(Ns):: Eigenv
    real(kind=Rkind),dimension(Ns, Ns):: bmatrixc, bmatrixa, qzmatrix, newqzmatrix, qxmatrix, newqxmatrix, b0matrixa, newb0matrixa, b0matrixc, newb0matrixc, x2matrix, newx2matrix, nmatrix, newnmatrix

    real(kind=Rkind),dimension(Dim):: vec_mid1, vec_mid2
    real(kind=Rkind),dimension(Dim, Dim):: Hn

    integer:: cevec=1                          
    integer:: ierr                                       
    real(kind=Rkind),dimension(Dim):: eigvalue, fv1, fv2
    real(kind=Rkind),dimension(Dim, Dim):: eigvector

    do i=1, Ns
        Eigenv(i)=Lamda*Eigenv(i)
    end do

    Hn=0.0_Rkind

!----------loc term-----------
    do i=1, Nb
        do j=1, Ns
            Hn(Ns*(i-1)+j, Ns*(i-1)+j)=Hn(Ns*(i-1)+j, Ns*(i-1)+j)+Eigenv(j)
        end do
    end do

!----------coupling term-----------
    bmatrixan1=0.0_Rkind
    bmatrixcn1=0.0_Rkind
    do i=1, Nb-1
        tmp=i
        tmp=sqrt(tmp)
        bmatrixan1(i, i+1)=bmatrixan1(i, i+1)+tmp
        bmatrixcn1(i+1, i)=bmatrixcn1(i+1, i)+tmp
    end do

    do i=1, Nb-1
        do j=1, Ns
            do k=1, Ns
            Hn(Ns*i+j, Ns*(i-1)+k)=Hn(Ns*i+j, Ns*(i-1)+k)+tn*bmatrixa(j, k)*bmatrixcn1(i+1, i)
            Hn(Ns*(i-1)+k, Ns*i+j)=Hn(Ns*(i-1)+k, Ns*i+j)+tn*bmatrixa(j, k)*bmatrixcn1(i+1, i)
            end do
        end do
    end do

!----------bath term-----------
    do i=1, Nb
        do j=1, Ns
            Hn(Ns*(i-1)+j, Ns*(i-1)+j)=Hn(Ns*(i-1)+j, Ns*(i-1)+j)+epin*(i-1)
        end do
    end do

!-------------- diagonalize Hamilt matrix ------------------
!---- Call RS subroutine --------    
!--using slatec in 1976: slower, but independent of compiler ---- 
!          call  rs (Dim, Dim, Hn, eigvalue, &
!                          cevec, eigvector, fv1, fv2, ierr)
!--using Lapack: faster, but must with intel + Lapack ----
          call  rs_cpu (Dim, Dim, Hn, eigvalue, &
                                  cevec, eigvector, fv1, fv2, ierr)
!---------------------
    if (ierr/=0) then
        print *, 'error'
        stop
    endif
!----------------------------------------------------

    do i=1, Ns
        Eigenv(i)=eigvalue(i)
    end do

    tmp=Eigenv(1)
    do i=1, Ns
        Eigenv(i)=Eigenv(i)-tmp
    end do

    do i=1, Ns

        vec_mid1=0.0_Rkind
        do j=1, Nb-1
            tmp=j
            tmp=sqrt(tmp)
            do k=1, Ns
            vec_mid1(Ns*(j-1)+k)=tmp*eigvector(Ns*j+k, i)
            end do
        end do
        
        do j=1, Ns

            sum1=0.0_Rkind
            do k=1, Dim
                sum1=sum1+eigvector(k, j)*vec_mid1(k)
            end do

            bmatrixa(j, i)=sum1
        end do
    end do


    do i=1, Ns
        do j=1, Ns
            bmatrixc(i, j)=bmatrixa(j, i)
        end do
    end do

!---- Write qx matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*qxmatrix(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newqxmatrix(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            qxmatrix(i, j)=newqxmatrix(i, j)
        end do
    end do

!---- Write qz matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*qzmatrix(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newqzmatrix(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            qzmatrix(i, j)=newqzmatrix(i, j)
        end do
    end do

!    do i=1, Ns
!        do j=1, Ns
!            if(abs(qzmatrix(i, j)-qzmatrix(j, i))>=0.00000001) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

!---- Write x^2 matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*x2matrix(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newx2matrix(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            x2matrix(i, j)=newx2matrix(i, j)
        end do
    end do

!    do i=1, Ns
!        do j=1, Ns
!            if(abs(x2matrix(i, j)-x2matrix(j, i))>=0.00000001) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

!---- Write <a+a> matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*nmatrix(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns

            sum2=0.0_Rkind
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newnmatrix(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            nmatrix(i, j)=newnmatrix(i, j)
        end do
    end do

!    do i=1, Ns
!        do j=1, Ns
!            if(abs(nmatrix(i, j)-nmatrix(j, i))>=0.00000001) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

!---- Write b0 matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*b0matrixa(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns
            sum2=0
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newb0matrixa(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            b0matrixa(i, j)=newb0matrixa(i, j)
        end do
    end do

!---- Write b0+ matrix --------
    do i=1, Ns

        vec_mid2=0.0_Rkind
        do j=1, Nb
            do k=1, Ns
                tmp=0.0_Rkind
                do l=1, Ns
                    tmp=tmp+eigvector(Ns*(j-1)+l ,i)*b0matrixc(k ,l)
                end do
                vec_mid2(Ns*(j-1)+k)=tmp
            end do
        end do

        do j=1, Ns
            sum2=0.0_Rkind
            do k=1, Dim
                sum2=sum2+eigvector(k, j)*vec_mid2(k)
            end do
            newb0matrixc(j, i)=sum2
        end do
    end do

    do i=1, Ns
        do j=1, Ns
            b0matrixc(i, j)=newb0matrixc(i, j)
        end do
    end do

!    do i=1, Ns
!        do j=1, Ns
!            if(abs(b0matrixa(i, j)-b0matrixc(j, i))>=0.000000001) then
!                print *, "e"
!                pause
!            end if
!        end do
!    end do

    end subroutine diag_Hn


!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: RS_CPU
! NOTICE :
! TYPE   : subroutine
! PURPOSE: This code use the same interface as RS() (in slatec), but
!          actually called the Lapack code "dsyev()" to diagonalize
!          the real symmetric matrix.
! I/O    : 
! VERSION: 26-Dec-2012
! AUTHOR : Q. S. Wu
! COMMENT: I got it from Quan-Shen Wu.
!========+=========+=========+=========+=========+=========+=========+=$
      SUBROUTINE RS_cpu (NM, N, A, W, MATZ, Z, FV1, FV2, IERR)

         use prec
         implicit none
    
         !>> inout variables 
         integer,intent(IN):: N, NM, MATZ
         integer,intent(OUT):: IERR
         real(kind=rkind),dimension(NM,NM),intent(IN):: A
         real(kind=rkind),dimension(NM,N),intent(OUT)::Z
         real(kind=rkind),dimension(N),intent(OUT)::W
         real(kind=rkind),dimension(N),intent(IN)::FV1,FV2

         !>> local variables
         integer :: i
         real(kind=rkind),allocatable :: ham(:, :)

         allocate(ham(N, N))

         ham= A(1:N, 1:N)

         call eigensystem_r('V', 'U', N, ham, W)

         Z(1:N, 1:N)= ham

         ierr=0

         return
      end subroutine rs_cpu
!========+=========+=========+=========+=========+=========+=========+=$



!========+=========+=========+=========+=========+=========+=========+=$
! PROGRAM: eigensystem_r
! NOTICE :
! TYPE   : subroutine
! PURPOSE: This code encapsulate the Lapack code "dsyev()" for diagonalizing
!          real symmetric matrix and prudce eigen values and eigen vectors.
! I/O    : 
! VERSION: 26-Dec-2012
! AUTHOR : Q. S. Wu
! COMMENT: I got it from Quan-Shen Wu.
!========+=========+=========+=========+=========+=========+=========+=$
! real version on cpu
! a subroutine to calculate eigenvector and eigenvalue
  subroutine eigensystem_r(JOBZ,UPLO,N,A,W)

!    use param, only : Rkind
     implicit none

     integer, parameter :: Rkind=8
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1, intent(in) :: JOBZ

!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
     character*1, intent(in) :: UPLO

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

     real(Rkind),intent(inout) :: A(N,N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

     real(Rkind), intent(inout) :: W(N)
    
     integer :: info

     integer :: lwork

     real(Rkind),allocatable ::  rwork(:)

     real(Rkind),allocatable :: work(:)

     lwork= 16*N
     allocate(rwork(lwork))
     allocate( work(lwork))


     W= 0.0d0
     info= 0
     call dsyev( JOBZ, UPLO, N, A, N,  &
              W, work, lwork, rwork, info )

     return
  end subroutine
!========+=========+=========+=========+=========+=========+=========+=$




