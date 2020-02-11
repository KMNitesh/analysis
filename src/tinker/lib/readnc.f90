! Created by ${USER_NAME} on 4/11/18.

subroutine readnetcdf(inetcdf)
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
    integer(kind = 8) inetcdf
    integer  ret;
    call netcdf_read_next(inetcdf,xbox,ybox,zbox,alpha,beta,gamma,x,y,z,ret)
    if (ret /= 1) then
        abort = .true.
        n = 0
        return
    end if
    abort = .false.

    call lattice

    return
end subroutine readnetcdf

subroutine opennetcdf(netcdffile, inetcdf)
    character*120 netcdffile
    integer(kind = 8) inetcdf
    integer leng,trimtext
    integer ret

    leng = trimtext(netcdffile)
    call netcdf_open(netcdffile, leng, inetcdf,ret)
    if (ret /= 0) then
        write(iout,*) 'error open file '//netcdffile
        stop
    end if

end subroutine opennetcdf

subroutine closenetcdf(inetcdf)
    integer(kind = 8) inetcdf
    call netcdf_close(inetcdf)

end subroutine closenetcdf

subroutine rewindnetcdf(inetcdf)
    integer(kind = 8) inetcdf
    call netcdf_rewind(inetcdf)

end subroutine rewindnetcdf

subroutine getnetcdf(netcdffile, usenetcdf)
    use inform
    use iounit
    use output
    implicit none
    integer(kind = 8) inetcdf
    logical exist
    character*120 netcdffile
    logical usenetcdf

    usenetcdf = .false.

    call nextarg (netcdffile, exist)
    if (exist) then
        if (netcdffile(1 : 1) /= 'Y' .and. netcdffile(1 : 1) /= 'y') then
             return
        end if
         usenetcdf = .true.
        call nextarg (netcdffile, exist)

        if (.not. exist) then
            goto 100
        end if
        call trimhead(netcdffile)

        return
    end if
    write (iout, 9)
    9 format (/, ' Use Amber NetCDF File : ', $)
    read (input, 20) netcdffile
    if (netcdffile(1 : 1) == 'Y' .or. netcdffile(1 : 1) == 'y') then
        usenetcdf = .true.
        100       do while (.not. exist)
            write (iout, 10)
            10     format (/, ' Enter NetCDF Coordinate File Name :  ', $)
            read (input, 20)  netcdffile
            20    format (a120)
             call trimhead(netcdffile)
            inquire (file = netcdffile, exist = exist)
        end do
    end if

    return
end subroutine getnetcdf
