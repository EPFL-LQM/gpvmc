        SUBROUTINE zdotu_sub(N,x,incx,y,incy,dotu)
            INTEGER N, incx, incy
            DOUBLE COMPLEX x(*), y(*), dotu, zdotu
            external zdotu
            dotu=zdotu(N,x,incx,y,incy)
            return
        end

