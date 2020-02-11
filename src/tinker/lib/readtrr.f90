! Created by ${USER_NAME} on 4/11/18.

subroutine readtrr(itrr)
    use sizes
    use atomid
    use atoms
    use bound
    use boxes
    use couple
    use files
    use inform
    use iounit
    use titles
    integer(kind = 8) itrr
    integer  ret;
    call trr_read_next(itrr, xbox, ybox, zbox, alpha, beta, gamma, x, y, z, ret)
    if (ret /= 1) then
        abort = .true.
        n = 0
        return
    end if
    abort = .false.

    call lattice

    return
end subroutine readtrr

subroutine opentrr(trrfilename, itrr)
    character*120 trrfilename
    integer(kind = 8) itrr
    integer leng,trimtext

    leng = trimtext(trrfilename)
    call trr_open(trrfilename, leng, itrr)

end subroutine opentrr

subroutine closetrr(itrr)
    integer(kind = 8) itrr
    call trr_close(itrr)

end subroutine closetrr

subroutine gettrr(trrfile, usetrr)
    use inform
    use iounit
    use output
    implicit none
    integer(kind = 8) itrr
    logical exist
    character*120 trrfile
    logical usetrr

    usetrr = .false.
    call nextarg (trrfile, exist)
    if (exist) then
        if (trrfile(1 : 1) /= 'Y' .and. trrfile(1 : 1) /= 'y') then
            return
        end if
        usetrr = .true.
        call nextarg (trrfile, exist)

        if (.not. exist) then
            goto 100
        end if
        call trimhead(trrfile)

        return
    end if
    write (iout, 9)
    9 format (/, ' Use Gromacs Trr File : ', $)
    read (input, 20) trrfile
    if (trrfile(1 : 1) == 'Y' .or. trrfile(1 : 1) == 'y') then
        usetrr = .true.
        100       do while (.not. exist)
            write (iout, 10)
            10     format (/, ' Enter Trr Coordinate File Name :  ', $)
            read (input, 20)  trrfile
            20    format (a120)
             call trimhead(trrfile)
            inquire (file = trrfile, exist = exist)
        end do
    end if

    return
end subroutine gettrr
