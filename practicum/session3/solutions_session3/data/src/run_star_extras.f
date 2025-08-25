! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def
!      use mlt_def
      
      implicit none

      real(dp) :: rmin
      real(dp) :: r_div_rcore
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         
         s% Abundance_pgstar_decorator => Abundance_pgstar_decorator
         s% HR_pgstar_decorator => HR_pgstar_decorator  
         s% History_Panels1_pgstar_decorator => History_Panels1_pgstar_decorator
         s% History_Panels2_pgstar_decorator => History_Panels2_pgstar_decorator

      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
	 !extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
         !At startup, give a value to the radius rzero
         rmin = s% r(1)
      end subroutine extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      ! Make a function for the pgplot decorator
      subroutine Abundance_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not dp
         real,intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num 
         real :: xcenter,ycenter,dx,dy,a
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s

         ! I think this part is assigning id the star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dx=(xmax-xmin)
         dy=(ymax-ymin)

         xcenter=xmin+dx*0.5
         ycenter=ymin+dy*0.5

         call pgsci(clr_CadetBlue)

         do i=1,4
            a=(i/10.0)
            call pgline(5, (/xcenter-a*dx,xcenter-a*dx,xcenter+a*dx,xcenter+a*dx,xcenter-a*dx/),&
                           (/ycenter-a*dy,ycenter+a*dy,ycenter+a*dy,ycenter-a*dy,ycenter-a*dy/))
         end do

         call pgsci(clr_Lilac)

         call pgptxt(xcenter,ycenter, 0.0, 1.0, 'Some added text on this plot')

      end subroutine Abundance_pgstar_decorator

      ! Make a function for the pgplot decorator for the HRD 
      !     This is basically putting the lines of constant radii
      subroutine HR_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not dp
         real,intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num
         real :: xcenter,ycenter,dx,dy,a
         real :: x1, x2, y1, y2, R, xtext, ytext
         real :: sigma_SB_SI, RSun_SI, LSun_SI
         Character(50) Rstring, Rtmp
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s

         ! I think this part is assigning id the star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Constants and values
         sigma_SB_SI = 5.670373e-8
         RSun_SI = 6.957e8
         LSun_SI = 3.846e26
 
         ! log(Teff) of min and max of the lines
         x1 = 0
         x2 = 6

         ! Calculate mid of plot
         dx=(xmax-xmin)
         dy=(ymax-ymin)
         xcenter=xmin+dx*0.5
         ycenter=ymin+dy*0.5

         call pgsci(clr_Gray)

         do i=-1,4,1
            ! Get the radius of the line
            R = (10**float(i))*RSun_SI
            ! Compute the functions for log L = k* log Teff + m
            y1 = log10(10**(log10(4*pi*(R**2) * sigma_SB_SI) + 4*x1)/LSun_SI)
            y2 = log10(10**(log10(4*pi*(R**2) * sigma_SB_SI) + 4*x2)/LSun_SI)

            ! Draw the line
            call pgline(2, (/x1,x2/),(/y1,y2/))
            
            ! Mark the constant radius close to the line
            ytext = ymin+0.05*dy
            xtext = (ytext - log10(4*pi*(R**2)*sigma_SB_SI/LSun_SI))/4 - 0.03*dx
            if (xtext < xmax) then
               xtext = xmax - 0.03*dx
               ytext = log10(10**(log10(4*pi*(R**2) * sigma_SB_SI) + 4*xtext)/LSun_SI) - 0.02*dy
            end if             

            write(Rtmp,'(f10.1)') 10**float(i)
            Rstring = trim(Rtmp) // trim(' R\d\(2281)')
            call pgptxt(xtext,ytext,0.0,1.0,Rstring)
         end do


      end subroutine HR_pgstar_decorator

      ! Make a function for the pgplot decorator for the logR - t diagram
      !     This is putting labels for case A, B and C
      subroutine History_Panels1_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not dp
         real,intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num
         real :: xcenter,ycenter,dx,dy,a
         real :: xA, yA, xB, yB, xC, yC
         real :: sigma_SB_SI, RSun_SI, LSun_SI
         Character(50) Rstring, Rtmp
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s

         ! I think this part is assigning id the star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dx = (xmax - xmin)
         dy = (ymax - ymin)

         ! Put case A
         if (s% center_h1 > 0.01) then
            xA = xmax
            yA = s% log_surface_radius
         end if
         call pgsci(clr_RoyalPurple)
         call pgline(2,(/xA-0.2*dx,xA/),(/yA,yA/))
         call pgptxt(xA-0.02*dx,yA+0.02*dy,0.0,1.0,'case A')

         ! Put case B
         if ((s% center_h1 < 0.01) .and. (s% center_he4 > 0.95)) then
            xB = xmax
            yB = s% log_surface_radius
         end if
         if (s% center_h1 < 0.01) then
            call pgsci(clr_LightSkyGreen)
            call pgline(2,(/xB-0.2*dx,xB/),(/yB,yB/))
            call pgptxt(xB-0.02*dx,yB+0.02*dy,0.0,1.0,'case B') 
         end if

         ! Put case C
         if ((s% center_h1 < 0.01) .and. (s% center_he4 < 0.95)) then
            xC = xmax
            yC = s% log_surface_radius
            call pgsci(clr_IndianRed)
            call pgline(2,(/xC-0.2*dx,xC/),(/yC,yC/))
            call pgptxt(xC-0.02*dx,yC+0.02*dy,0.0,1.0,'case C')           
         end if

      end subroutine History_Panels1_pgstar_decorator

      ! Make a function for the pgplot decorator for the P-t diagram
      !     This is putting labels for the mass transfer cases A, B and C
      subroutine History_Panels2_pgstar_decorator(id, xmin, xmax, ymin, ymax, plot_num, ierr)
         integer, intent(in) :: id
         !Not dp
         real,intent(in) :: xmin, xmax, ymin, ymax
         integer, intent(in) :: plot_num
         real :: dx,dy, G, q, M2, rL, separation, period_days
         real :: xA, yA, xB, yB, xC, yC, xq, yq
         real :: sigma_SB_SI, RSun_SI, LSun_SI
         Character(50) qstring, qtmp
         integer, intent(out) :: ierr
         integer :: i
         type (star_info), pointer :: s

         ! I think this part is assigning id the star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dx = (xmax - xmin)
         dy = (ymax - ymin)

         G = 4*pi*pi/(365.25*365.25)     ! AU^3 MSun^-1 days^-2
         q = s% x_ctrl(1)
         M2 = q * s% star_mass
         rL = 0.49 * q**(2/3) / (0.69 * q**(2/3) + log(1 + q**(1/3)))
         separation = (10**s% log_surface_radius / 215) / rL
         period_days = sqrt(separation**3 * 4*pi*pi/(G*s% star_mass * M2))

         ! Mark the mass ratio assumed
         write(qtmp,'(f10.3)') q
         qstring = trim('Assuming q = ') // trim(qtmp)
         xq = xmax - 0.05*dx
         yq = ymax - 0.05*dy
         call pgptxt(xq,yq,0.0,1.0,qstring)

         ! Put case A
         if (s% center_h1 > 0.01) then
            xA = xmax
            yA = log10(period_days)
         end if
         call pgsci(clr_RoyalPurple)
         call pgline(2,(/xA-0.2*dx,xA/),(/yA,yA/))
         call pgptxt(xA-0.02*dx,yA+0.02*dy,0.0,1.0,'case A')

         ! Put case B
         if ((s% center_h1 < 0.01) .and. (s% center_he4 > 0.95)) then
            xB = xmax
            yB = log10(period_days)
         end if
         if (s% center_h1 < 0.01) then
            call pgsci(clr_LightSkyGreen)
            call pgline(2,(/xB-0.2*dx,xB/),(/yB,yB/))
            call pgptxt(xB-0.02*dx,yB+0.02*dy,0.0,1.0,'case B')
         end if

         ! Put case C
         if ((s% center_h1 < 0.01) .and. (s% center_he4 < 0.95)) then
            xC = xmax
            yC = log10(period_days)
            call pgsci(clr_IndianRed)
            call pgline(2,(/xC-0.2*dx,xC/),(/yC,yC/))
            call pgptxt(xC-0.02*dx,yC+0.02*dy,0.0,1.0,'case C')
         end if

      end subroutine History_Panels2_pgstar_decorator


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: log10_period_days, period_days, separation, G, rL, M2, q
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         names(1) = "period_days"
         G = 4*pi*pi/(365.25*365.25)     ! AU^3 MSun^-1 days^-2
         q = s% x_ctrl(1)
         M2 = q * s% star_mass
         rL = 0.49 * q**(2/3) / (0.69 * q**(2/3) + log(1 + q**(1/3)))
         separation = (10**s% log_surface_radius / 215) / rL
         period_days = sqrt(separation**3 * 4*pi*pi/(G*s% star_mass * M2))
         vals(1) = period_days

         names(2) = "log_period_days"
         vals(2) = log10(period_days)

      end subroutine data_for_extra_history_columns


      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
     
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

         ! we want the relaxed rmin, not the r at the first timestep
         rmin = min(rmin, s% r(1))
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info

      end module run_star_extras
      

