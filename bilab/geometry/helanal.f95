!*******************************************************************************
! Description:
!
!
! Reference:
!  Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics
! of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586.
!  
!  Sugeta H, Miyazawa T. General method for calculating helical parameters of
! polymer chains from bond lengths, bond angles, and internal-rotation angles.
! Biopolymers. 1967;5: 673–679. 
!
! Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!
!
!  Created:
!    11 Janurary 2016
!
!  Author:
!
!    Wei Cao
!
! Compilation: (DEBUG)
!  gfortran:
!     gfortran -g -cpp -DQUAD_PRECISION -c helanal.f95
!     gfortran -g -DQUAD_PRECISION -m64 -fbounds-check -Wall -fbacktrace -Wtabs\ 
!            -finit-real=nan \
!            test_f95.f95 helanal.o -o test_f95
!  f2py: manually generate python extension or see setup.py 
!    f2py-2.7 helanal.f95 -m _helanal -h helanal.pyf --overwrite-signature
!    f2py-2.7 -c helanal.pyf helanal.f95
!*******************************************************************************
!****************************************************************************
!                                                                           *
! PROGRAM TO CHARACTERISE THE GEOMETRIES OF ALPHA HELICES IN PROTEINS.      *
!                                                                           *
! AUTHORS: Sandeep Kumar and Prof. Manju Bansal,                            *
!          MBU, Indian Institute of Science, Bangalore 560012, India        *
!                                                                           *
!  e-mail address: mb@mbu.iisc.ernet.in                                     *
!                                                                           *
! Method outlined in the following paper:                                   *
!                                                                           *
! Kumar, S. and Bansal, M. (1996). Structural and Sequence Characteristics  * 
! of Long Alpha Helices in Globular Proteins. Biophysical J.,71, 1574-1586. *
!                                                                           *
! If any Bug is found, please report it to the authors.                     *
!                                                                           *
! SUMMARY OF THE ALGORITHM                                                  *
!                                                                           *
! Geometry of an alpha helix is characterised in terms of the angles between*
! local helix axes and the path traced by the local helix origins,          *
! calculated using the procedure of Sugeta and Miyazawa (1967), for every   *
! set of four contiguous C-alpha atoms, and sliding this window over the    *
! length of the helix in steps of one C-alpha atom.                         *
! Matrix M(I, J) contains the bending angles between local helix axes I & J *
! which are used to characterize the overall geometry of the  helix.        * 
! The local helix origins trace out the path described by the helix in 3-D  *
! space. These origins are reoriented in X-Y plane and the reoriented       *
! points are used to fit a circle and a line by least squares method.       *
! Unit twist and unit height of the alpha helix are also calculated.        *
!                                                                           *
! A maximum of 5000 helices, each with 100 or less number of residues can   *
! be analysed at one time.                                                  * 
!                                                                           *
! In order to analyse a larger number of helices, increase the value of 'i' *
! in the  following statement:                                              *
!                                                                           *
!      parameter (i=5000) in the main program.                              *
!                                                                           *
! In order to analyse helices with lengths greater than 100 residues,       *
! increase the dimensions of the appropriate variables in the dimension     *
! statements of the main program and subroutines.                           *
!                                                                           *
! This program has several options for input and is fully interactive.      *
!                                                                           *
! INPUT FILES :                                                             *
!                                                                           *
! Input files to this program can be different depending upon the options   *
! chosen. This program can directly read HELIX records in one or more PDB   *
! files. In order to analyse helices found in helix records of a PDB file,  *
! give the PDB file name as input to the program. In order to analyse the   *
! helices found in the HELIX records of more than one PDB  files, give the  *
! PDB file names sequentially when prompted by the program or input the     *
! name of the file containing the PDB files names in format (5x,a11) as the *
! input to the program.                                                     *
!                                                                           *
! HELIX records can be read from PDB files or can be read from a file       *
!      containing information about the helix start/end residues written in *
!      the same format as PDB HELIX records or in a different format, in    *
!      which case the format has to be keyed in when running the program.   *
! In the last case, files other than the PDB files or files containing      *
!      C-alpha coordinates in the PDB format, can also be given as input to *
!      the program, upon specifying their format.                           * 
!                                                                           *
!                                                                           *
! All PDB files and other input files (if any) should be in the same        *
! directory.                                                                *
!                                                                           *
! OUTPUT FILES :                                                            *
!                                                                           *
! For each  run of HELANAL on file(s) containing alpha helices with length  *
! greater than or equal to 9 residues, the following output files are       *
! created.                                                                  *
!                                                                           *
! RUN.ANS contains the questions and their answers during a run of HELANAL. *
!                                                                           *
! HELINFO.OUT file created only when the HELIX records in the PDB files are *
! used for information on helix start/end residues, contains these helix    *
! records.                                                                  *
!                                                                           *
! HELCA.OUT contains Coordinates of the C-alpha atoms constituting the      *
! helices.                                                                  *
!                                                                           *
! AXES.OUT contains the local helix axes fitted to 4 consecutive C-alpha    *
! atoms along with a matrix  M(I, J) whose elements are the angles between  *
! local helix axes I and J.                                                 *
!                                                                           *
! ANGLE.OUT contains the angle between successive local helix axes. It also *
! lists mean bending angle and maximum bending angle for each helix.        *
!                                                                           *
! ORIGIN.OUT contains the local helix origins for the helix along with the  *
! statistics obtained by fitting least square plane, circle and line to the *
! local helix origins.                                                      *
!                                                                           *
! NH.OUT contains unit height and unit twist for every turn of the helix as *
! well as average unit height and average unit twist for the whole helix.   *  
!                                                                           *
! ***.PRM contains summary of various parameters obtained by HELANAL.       *
!                                                                           *
! ***.TAB contains the parameters from ***.PRM file in a tabular form along *
! the overall geometry assignment.                                          *
!                                                                           *
!****************************************************************************
module helanal
  ! *-- Global --*
  implicit none
  ! *-- symbolic name: f2py could not recognize them if used in sub/func --*
  ! *-- f2py not recognize them --* 
!#ifdef QUAD_PRECISION
!  integer, parameter :: DP = 16 ! kind(1.0d0)
!#elif TEN_DIGIT_PRECISION
!  integer, parameter :: DP = selected_real_kind(10) !kind(1.0d0)
!#else
!  integer, parameter :: DP = 8 ! kind(1.0d0)
!#endif

  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
  ! *-- const parameters --*
  real(DP), parameter :: pi = 3.141592653589793238462643383279502884197
  ! *-- public: default is private --*
  public :: outer
  public :: cross
  public :: fit
  public :: local_helix
  !public :: fit_origins_lsq

contains

  !*****************************************************************************
  ! Description:
  !   print 2d matrix for a given points (array)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status: 
  !    OK
  ! Arguments:
  !   mat (array): Nxd elements, N samples and d dimensions
  !
  ! Returns:
  !   None
  !*****************************************************************************
  subroutine mat2d_print(mat)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:, :), intent(in):: mat
    integer :: nrows, i
    nrows = size(mat, dim=1)
    do i=1, nrows
       print *, mat(i, :)
    end do
  end subroutine mat2d_print

  !*****************************************************************************
  ! Description:
  !   Calculate the outer product of two vectors, a and b
  !   (Not used now)
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK 
  ! Arguments:
  !   a (array): a vector of N elements
  !   b (array): a vector of N elements
  !
  ! Returns:
  !    N x N outer product
  !*****************************************************************************
  function outer(a, b)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    ! description of input and output
    real(DP), dimension(:), intent(in) :: a, b
    real(DP), dimension(size(a), size(b)) :: outer
    !f2py real(DP) intent(hide), depend(a, b):: outer

    outer = spread(a, dim=2, ncopies=size(b)) * &
         spread(b, dim=1, ncopies=size(a))

  end function outer


  !*****************************************************************************
  ! Description:
  !    Calculate cross product of two vectors, a and b
  !
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK 
  !
  ! Arguments:
  !   a (array): a vector of 3 elements
  !   b (array): a vector of 3 elements
  !
  ! Returns:
  !   a unit vector of 3 elements, which is perpendicular
  !   to the ab-plane
  !*****************************************************************************
  function cross(a, b)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(3), intent(in) :: a, b
    real(DP), dimension(3) :: cross
    real(DP) :: length

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

    length = sqrt(sum(cross * cross))
    ! normalization
    cross = cross * 1.0 / length
  end function cross

  !*****************************************************************************
  ! Description:
  !   For a given four continuous 3d points, find direction along the helix 
  !   axis, origins, twist angle, and height of the turn of a helix
  !  
  ! Standard format:
  !    f90 and later
  !
  ! Status:
  !    OK 
  !
  ! Arguments:
  !   p1 (array): a vector of 1x3 array
  !   p2 (array): a vector of 1x3 array
  !   p3 (array): a vector of 1x3 array
  !   p4 (array): a vector of 1x3 array
  !
  ! Returns:
  !   direct
  !   origin
  !   twist
  !   height
  !*****************************************************************************
  subroutine local_helix(p1, p2, p3, p4, direct, origin, twist, height)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:), intent(in):: p1, p2, p3, p4
    real(DP), dimension(:), intent(out):: direct
    real(DP), dimension(:, :), intent(out):: origin
    real(DP), intent(out):: twist, height

    ! local variables
    real(DP), dimension(3) :: v12, v23, v34, dv13, dv24
    real(DP) :: length1, length2, r, cos_theta

    v12 = p2 - p1
    v23 = p3 - p2
    v34 = p4 - p3

    dv13 = v12 - v23
    dv24 = v23 - v34

    ! normal vector
    direct = cross(dv13, dv24)

    ! angle
    length1 = sqrt(dot_product(dv13, dv13))
    length2 = sqrt(dot_product(dv24, dv24))
    cos_theta = dot_product(dv13, dv24) / (length1 * length2)

    ! twist: in degree
    twist = acos(cos_theta) * 180.0 / pi

    ! radius of local helix cylinder
    r = sqrt(length1 * length2) / 2.0 *(1 - cos_theta)

    ! height of local helix cylinder
    height = dot_product(v23, direct)

    ! unit vectors
    dv13 = dv13 / norm2(dv13)
    dv24 = dv24 / norm2(dv24)

    ! origins
    origin(1, :) = p2 - dv13 * r
    origin(2, :) = p3 - dv24 * r

    ! rotation vector
    ! rot_vectors(1, :) = dv13
    ! rot_vectors(2, :) = dv24

  end subroutine local_helix

  !*****************************************************************************
  ! Description:
  !   Fit least-squared plane to local helix origins (Not finished )
  !   (not used now)
  ! Standard format:
  !   f90 and later
  !
  ! Arguments:
  !   origins (array): a vector of 1x3 array
  !
  ! Returns:
  !   radc
  !   rmsdc
  !   rmsdl
  !   r2
  !*****************************************************************************
  subroutine fit_origins_lsq(origins, radc, rmsdc, rmsdl, r2)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:,:), intent(in):: origins
    real(DP), intent(out):: radc, rmsdc, rmsdl, r2

    ! *-- local variables --*
    integer :: nrows, ncols
    real(DP), dimension(size(origins, 2)):: centroid ! dim=2: get columns
    real(DP), dimension(size(origins, 1), size(origins, 2)):: M

    nrows = size(origins, 1)
    ncols = size(origins, 2)
    centroid = sum(origins, 2) / ncols

    ! centralize
    M = origins - spread(centroid, dim=1, ncopies=nrows)
    ! *-- SVD --*
    
  end subroutine fit_origins_lsq


  !*****************************************************************************
  ! Description:
  !   For a set of CA coordinates, find bending angles along the helix axis
  !   using a sliding window of 9 points.
  !
  ! Standard format:
  !   f90 and later
  !
  ! Status:
  !   partially completed. fit_origins_lsq will be completed in future.
  ! 
  ! Arguments:
  !   points (array): Nx3 matrix
  !
  ! Returns:
  !   bending_angles (array):
  !   r2 (real(DP)) : estimated radius
  !*****************************************************************************
  subroutine fit(points, bending_angles)
    integer, parameter :: SP = kind(1.0)
    integer, parameter :: DP = kind(1.0d0)
    real(DP), dimension(:, :), intent(in):: points
    real(DP), dimension(shape(points,1)), intent(out):: bending_angles
    !f2py real(DP), intent(out) bending_angles
    ! real(DP), intent(out):: radc, rmsdc, rmsdl
    ! real(DP), intent(out):: r2

    ! *-- local variables --*
    integer :: nrows, ncols, i, dsize
    integer, dimension(2):: s
    !real(DP), dimension(3, 3):: rotMat
    real(DP), dimension(3) :: direct
    real(DP), dimension(2, 3):: origin
    real(DP) :: twist, height, angle, tmp

    ! *-- dynamic local variables --*
    real(DP), dimension(:, :), allocatable :: directions, origins
    real(DP), dimension(:), allocatable :: twists, heights

    s = shape(points)
    nrows = s(1)
    ncols = s(2)
    dsize = ncols - 3

    ! *-- allocate memory --*
    allocate(directions(dsize, 3))
    allocate(origins(dsize, 3))
    allocate(twists(dsize))
    allocate(heights(dsize))

    ! *-- allocate(bending_angles(dsize - 3)) --*
    if (nrows.ne.3) then 
      stop 'Wrong shape of the input points'
    end if
    do i=1, dsize
       ! *-- slide along the helix spiral curve --*
       call local_helix(points(i,:), points(i+1,:), &
                        points(i+2, :), points(i+3,:), &
                        direct, origin, twist, height)
       ! *-- save parameters of local helix --*
       directions(i, :) = direct
       origins(i, :) = origin(1, :)
       origins(i+1, :) = origin(2, :)
       twists(i) = twist
       heights(i) = height
    end do
    ! *-- bending angle --*
    do i=1, dsize-3
       angle = dot_product(directions(i,:), directions(i+3, :))
       if (abs(angle - 1.0).le.1.0E-06) angle = 1.0
       ! *-- angle in degree --*
       bending_angles(i) = acos(angle) * pi
    end do

    ! *-- fit to circle --*
    ! reorient origins move to python
    ! *-- deallocate --*
    deallocate(directions)
    deallocate(origins)
    deallocate(twists)
    deallocate(heights)

  end subroutine fit

end module helanal
