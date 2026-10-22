! EXACT v2.2 (vectorized): Lambert solver, elliptic single-revolution transfers.
!   1. guess: z from a precomputed 64 x 32 lookup table (Taylor step near z = 0)
!   2. contour: one precomputed thin rectangle around the root, upper half only (symmetry),
!      8 Gauss-Legendre points on the top edge + 1 point per side; root = ratio of two weighted sums
!   3. one Newton correction; safety retry on the right rectangle if needed
! Same code as contour_study/gauss_legendre_vect.F90. Compile with:
!   gfortran -O3 -march=native -fno-trapping-math ... exact_vect.F90 ...
module rect_mod_fast
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: sqrt2 = 1.4142135623730950488016887_dp
    real(dp), parameter :: c_tp = 0.4714045207910316829338962_dp
    real(dp), parameter :: z_max = 4.0_dp * pi**2
    real(dp), parameter :: pio2_hi = 1.5707963267948966_dp, pio2_lo = 6.123233995736766e-17_dp   ! pi/2 = hi + lo
    logical :: initialized = .false.
    real(dp) :: nan_val

    integer, parameter :: N_tile = 40, N_v = 1, N_h = 8, NP = 2 * N_v + N_h
    real(dp), parameter :: h_tile = 1.0_dp, eps2 = 1.0e-7_dp
    integer, parameter :: j_lo = -40, j_hi = 1, NJ = j_hi - j_lo + 1
    real(dp), parameter :: z_thr = 2.0_dp
    ! Guess table: z on a grid of x = sqrt(y_end/r12) and v = sqrt(tp/dt), one table per sign of A
    ! (z depends only on tau = A/r12 and the normalized time; x stretches the r1 ~ r2 corner, v = 1 is z = 0)
    integer, parameter :: Nx_g = 64, Nv_g = 32
    real(dp) :: zg_tab(0:Nx_g-1, 0:Nv_g-1, 2)
    ! One-sqrt shortcut for sqrt(y) at contour points where y is nearly real and positive:
    ! if wr > r_1sqrt*|wi|, |y| = wr*sqrt(1 + (wi/wr)^2) is taken as wr (relative error < 1e-12)
    ! (serial build only; compile with -DRECT_VECTOR for the vectorized build, which uses the
    !  branch-free two-sqrt form so the contour loops vectorize across the block)
    real(dp), parameter :: r_1sqrt = 1.0e6_dp
    ! Number of problems solved side by side in solve_block
    integer, parameter :: NBLK = 32
    ! Refuse r1 ~ r2 (transfer angle within ~0.06 deg of 0 or 360 deg), the problem's singularity,
    ! with the same rule as ivLam: sqrt(2)/2 - |tau| < tau_gap, tau = A/r12. Since
    ! sqrt(2)/2 - |tau| = y_end/(sqrt(2)*r12), the test is y_end < sqrt(2)*tau_gap*r12.
    real(dp), parameter :: tau_gap = 1.0e-7_dp
    integer, parameter :: id_end0 = N_tile + 1, id_end4 = N_tile + 1 + NJ, N_all = N_tile + 1 + 2 * NJ

    real(dp), dimension(NP, 0:N_all-1) :: tK_re, tK_im, tq_re, tq_im
    real(dp), dimension(NP, 0:N_all-1, 2) :: tu_re, tu_im
    real(dp), dimension(NP, 0:N_all-1) :: tW_re, tW_im, tWz_re, tWz_im
    real(dp), dimension(0:N_all-1) :: tile_c, tile_a, tile_b

    integer, parameter :: M_sc = 1024
    real(dp), parameter :: h_sc = pi / real(M_sc, dp), inv_h_sc = real(M_sc, dp) / pi
    real(dp) :: tab_sin0(0:M_sc+1), tab_cos0(0:M_sc+1)

contains

    pure subroutine sincos_tab(a, sn, cs)
        real(dp), intent(in)  :: a
        real(dp), intent(out) :: sn, cs
        integer :: k
        real(dp) :: d, d2, sd, cd
        k = nint(a * inv_h_sc)
        d = a - real(k, dp) * h_sc
        d2 = d * d
        sd = d * (1.0_dp - d2 * (1.0_dp/6.0_dp - d2 * (1.0_dp/120.0_dp)))
        cd = 1.0_dp - d2 * (0.5_dp - d2 * (1.0_dp/24.0_dp - d2 * (1.0_dp/720.0_dp)))
        sn = tab_sin0(k) * cd + tab_cos0(k) * sd
        cs = tab_cos0(k) * cd - tab_sin0(k) * sd
    end subroutine sincos_tab

    subroutine gauss_legendre(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        integer :: i, k, it
        real(dp) :: t, p0, p1, p2, dp_n, dt

        do i = 1, (n + 1) / 2
            t = cos(pi * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
            do it = 1, 100
                p0 = 1.0_dp
                p1 = t
                do k = 2, n
                    p2 = ((2*k - 1) * t * p1 - (k - 1) * p0) / real(k, dp)
                    p0 = p1
                    p1 = p2
                end do
                dp_n = real(n, dp) * (t * p1 - p0) / (t * t - 1.0_dp)
                dt = p1 / dp_n
                t = t - dt
                if (abs(dt) < 1.0e-15_dp) exit
            end do
            x(i) = -t
            x(n + 1 - i) = t
            w(i) = 2.0_dp / ((1.0_dp - t * t) * dp_n * dp_n)
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre

    subroutine init_cache()
        integer :: j, k, m
        real(dp) :: zero, L

        do k = 0, N_tile
            call build_tile(k, real(k, dp) * z_max / N_tile, max(real(k, dp) * z_max / N_tile - h_tile, 0.0_dp), &
                            min(real(k, dp) * z_max / N_tile + h_tile, z_max), eps2)
        end do
        do j = j_lo, j_hi
            L = 2.0_dp**j
            call build_tile(id_end0 + j - j_lo, 1.5_dp * L, 0.25_dp * L, 3.0_dp * L, min(eps2, 0.1_dp * L))
            call build_tile(id_end4 + j - j_lo, z_max - 1.5_dp * L, z_max - 3.0_dp * L, z_max - 0.25_dp * L, &
                            min(eps2, 0.1_dp * L))
        end do

        call build_guess()

        do m = 0, M_sc + 1
            tab_sin0(m) = sin(real(m, dp) * h_sc)
            tab_cos0(m) = cos(real(m, dp) * h_sc)
        end do

        zero = 0.0_dp
        nan_val = zero / zero
        initialized = .true.
    end subroutine init_cache

    ! Normalized time of flight T(z) = dt*sqrt(mu)/r12^1.5 for y_end/r12 = yeh, tau = A/r12
    pure real(dp) function t_norm(z, yeh, tau)
        real(dp), intent(in) :: z, yeh, tau
        real(dp) :: s, y, C, Sv
        s = sqrt(z)
        if (tau > 0.0_dp) then
            y = yeh + 2.0_dp * sqrt2 * tau * sin(0.25_dp * s)**2
        else
            y = yeh - 2.0_dp * sqrt2 * tau * cos(0.25_dp * s)**2
        end if
        C = 2.0_dp * sin(0.5_dp * s)**2 / z
        if (z < 0.5_dp) then
            Sv = 1.0_dp/6.0_dp + z*(-1.0_dp/120.0_dp + z*(1.0_dp/5040.0_dp + z*(-1.0_dp/362880.0_dp &
               + z*(1.0_dp/39916800.0_dp + z*(-1.0_dp/6227020800.0_dp)))))
        else
            Sv = (s - sin(s)) / (s * z)
        end if
        t_norm = (y / C)**1.5_dp * Sv + tau * sqrt(y)
    end function t_norm

    ! Guess table by bisection on T(z) = T0/v^2 at each node (init only)
    subroutine build_guess()
        integer :: i, j, k, is
        real(dp) :: yeh, tau, y0, T0, Tt, lo, hi, zm
        do is = 1, 2
            do i = 0, Nx_g - 1
                yeh = max(1.0e-12_dp, (real(i, dp) / (Nx_g - 1))**2)
                tau = (1.0_dp - yeh) / sqrt2
                if (is == 2) tau = -tau
                y0 = yeh
                if (is == 2) y0 = yeh + 2.0_dp * sqrt2 * abs(tau)
                T0 = c_tp * y0 * sqrt(y0) + tau * sqrt(y0)
                zg_tab(i, 0, is) = z_max
                zg_tab(i, Nv_g - 1, is) = 0.0_dp
                do j = 1, Nv_g - 2
                    Tt = T0 / (real(j, dp) / (Nv_g - 1))**2
                    lo = 0.0_dp
                    hi = z_max
                    do k = 1, 60
                        zm = 0.5_dp * (lo + hi)
                        if (t_norm(zm, yeh, tau) > Tt) then
                            hi = zm
                        else
                            lo = zm
                        end if
                    end do
                    zg_tab(i, j, is) = 0.5_dp * (lo + hi)
                end do
            end do
        end do
    end subroutine build_guess

    ! sin(x) for |x| <= pi/2 by its Taylor series to x^21 (truncation < 4e-19), branch-free
    pure elemental real(dp) function sin_poly(x)
        real(dp), intent(in) :: x
        real(dp) :: x2
        x2 = x * x
        sin_poly = x * (1.0_dp + x2*(-1.0_dp/6.0_dp + x2*(1.0_dp/120.0_dp + x2*(-1.0_dp/5040.0_dp &
                 + x2*(1.0_dp/362880.0_dp + x2*(-1.0_dp/39916800.0_dp + x2*(1.0_dp/6227020800.0_dp &
                 + x2*(-1.0_dp/1307674368000.0_dp + x2*(1.0_dp/355687428096000.0_dp &
                 + x2*(-1.0_dp/121645100408832000.0_dp + x2*(1.0_dp/51090942171709440000.0_dp)))))))))))
    end function sin_poly

    ! Ratio of the top-edge sums sum c_i Im(z_i f_i), sum c_i Im(f_i), f = 1/Phi = q3/G, for one node set
    pure subroutine top_sums(n, u_re, u_im, K_re, K_im, q_re, q_im, P_re, P_im, Pz_re, Pz_im, &
                             y_end, s2absA, A, mu_t, num, den)
        integer, intent(in) :: n
        real(dp), intent(in) :: u_re(n), u_im(n), K_re(n), K_im(n), q_re(n), q_im(n)
        real(dp), intent(in) :: P_re(n), P_im(n), Pz_re(n), Pz_im(n), y_end, s2absA, A, mu_t
        real(dp), intent(out) :: num, den
        integer :: j
        real(dp) :: wr, wi, qr, qi, tt, ar, ai, br, bi, cr, ci, dr, di, inv_d
        num = 0.0_dp
        den = 0.0_dp
        do j = 1, n
            wr = y_end + s2absA * u_re(j)
            wi = s2absA * u_im(j)
            qr = sqrt(wr*wr + wi*wi) + wr                ! branch-free so the loop vectorizes
            tt = sqrt(2.0_dp * qr)
            qi = wi
            ar = wr*qr - wi*qi
            ai = wr*qi + wi*qr
            br = ar*K_re(j) - ai*K_im(j)
            bi = ar*K_im(j) + ai*K_re(j)
            cr = A*qr - tt*mu_t
            ci = A*qi
            dr = br + cr*q_re(j) - ci*q_im(j)
            di = bi + cr*q_im(j) + ci*q_re(j)
            inv_d = tt / (dr*dr + di*di + 1.0e-300_dp)
            num = num + (Pz_im(j)*dr - Pz_re(j)*di) * inv_d
            den = den + (P_im(j)*dr - P_re(j)*di) * inv_d
        end do
    end subroutine top_sums

    subroutine build_tile(id, cen, a, b, eps_t)
        integer, intent(in) :: id
        real(dp), intent(in) :: cen, a, b, eps_t
        integer :: j, side, m, n_seg
        real(dp) :: xv(N_v), wv(N_v), xh(N_h), wh(N_h), xg(max(N_v, N_h)), wg(max(N_v, N_h))
        complex(dp) :: corner(4), zp, dz, w, sh, chc, K1, q3

        call gauss_legendre(N_v, xv, wv)
        call gauss_legendre(N_h, xh, wh)
        tile_c(id) = cen
        tile_a(id) = a
        tile_b(id) = b
        corner(1) = cmplx(b, 0.0_dp, dp)
        corner(2) = cmplx(b, eps_t, dp)
        corner(3) = cmplx(a, eps_t, dp)
        corner(4) = cmplx(a, 0.0_dp, dp)
        m = 0
        do side = 1, 3
            if (side == 2) then
                n_seg = N_h
                xg(1:N_h) = xh
                wg(1:N_h) = wh
            else
                n_seg = N_v
                xg(1:N_v) = xv
                wg(1:N_v) = wv
            end if
            do j = 1, n_seg
                m = m + 1
                zp = 0.5_dp * (corner(side) + corner(side+1)) &
                   + 0.5_dp * (corner(side+1) - corner(side)) * xg(j)
                dz = 0.5_dp * (corner(side+1) - corner(side)) * wg(j)
                w = sqrt(zp)
                sh = sin(0.5_dp * w)
                chc = cos(0.5_dp * w)
                K1 = (w - 2.0_dp * sh * chc) / (2.0_dp * sqrt2 * w**3)
                q3 = (sh / w)**3
                tu_re(m, id, 1) = real(2.0_dp * sin(0.25_dp * w)**2, dp); tu_im(m, id, 1) = aimag(2.0_dp * sin(0.25_dp * w)**2)
                tu_re(m, id, 2) = real(2.0_dp * cos(0.25_dp * w)**2, dp); tu_im(m, id, 2) = aimag(2.0_dp * cos(0.25_dp * w)**2)
                tK_re(m, id) = real(K1, dp);   tK_im(m, id) = aimag(K1)
                tq_re(m, id) = real(q3, dp);   tq_im(m, id) = aimag(q3)
                tW_re(m, id) = real(dz, dp);   tW_im(m, id) = aimag(dz)
                tWz_re(m, id) = real((zp - cen) * dz, dp)
                tWz_im(m, id) = aimag((zp - cen) * dz)
            end do
        end do
    end subroutine build_tile

    subroutine rectsolution(R1, R2, dt, mu, H, V1, V2, Q)
        integer, intent(in) :: Q
        real(dp), intent(in) :: R1(3,Q), R2(3,Q), dt(Q), mu, H(3)
        real(dp), intent(out) :: V1(3,Q), V2(3,Q)
        integer :: k, nb
        real(dp) :: sqrt_mu, inv_sqrt_mu

        if (.not. initialized) call init_cache()

        sqrt_mu = sqrt(mu)
        inv_sqrt_mu = 1.0_dp / sqrt_mu

        ! Solve NBLK problems side by side: every step below is a loop over the block, so the
        ! long sqrt/division chains of independent problems overlap instead of running one at a time
        do k = 1, Q, NBLK
            nb = min(NBLK, Q - k + 1)
            call solve_block(R1(:,k:k+nb-1), R2(:,k:k+nb-1), dt(k:k+nb-1), sqrt_mu, inv_sqrt_mu, H, &
                             V1(:,k:k+nb-1), V2(:,k:k+nb-1), nb)
        end do
    end subroutine rectsolution

    ! Contour tile for a root estimate zz: scaled end tiles near 0 / 4*pi^2, regular tiles elsewhere
    pure integer function pick_tile(zz)
        real(dp), intent(in) :: zz
        if (zz < z_thr) then
            pick_tile = id_end0 + min(j_hi, max(j_lo, exponent(zz / 1.5_dp * sqrt2) - 1)) - j_lo
        else if (z_max - zz < z_thr) then
            pick_tile = id_end4 + min(j_hi, max(j_lo, exponent((z_max - zz) / 1.5_dp * sqrt2) - 1)) - j_lo
        else
            pick_tile = min(N_tile, max(0, nint(zz * N_tile / z_max)))
        end if
    end function pick_tile

    ! Single analytic polish step at g = zg; returns dz and y, sqrt(y) at g + dz
    pure subroutine polish(zg, ye, s2a, Ap, mt, dz, y, sy)
        real(dp), intent(in) :: zg, ye, s2a, Ap, mt
        real(dp), intent(out) :: dz, y, sy
        real(dp) :: sh, ch, sq, sh2, ch2, Cg, Sg, yg, syg, x3, Phi, dPhi
        real(dp) :: rinv, inv_sq, inv_sh2, inv_syg, inv_z, sqC, invC, dy
        ! Single analytic polish step z* = g - Phi(g)/Phi'(g),
        ! Phi = (y/C)^(3/2) S + A sqrt(y) - sqrt(mu) dt,
        ! Phi' = (y/C)^(3/2) [(2C^2 - 3S)/(4gC) + 3S^2/(4C)] + (A/8) [3S sqrt(y)/C + A sqrt(C)/sqrt(y)]
        ! One sincos at sqrt(g)/4; sin, cos at sqrt(g)/2 by double angle; sqrt(C) = sqrt(2) sin(sqrt(g)/2)/sqrt(g).
        sq = sqrt(zg)
        call sincos_tab(0.25_dp * sq, sh, ch)
        sh2 = 2.0_dp * sh * ch                           ! sin(sqrt(g)/2)
        ch2 = (ch - sh) * (ch + sh)                      ! cos(sqrt(g)/2)
        yg = ye + 2.0_dp * s2a * merge(sh * sh, ch * ch, Ap > 0.0_dp)   ! y(g) without cancellation
        yg = max(1.0e-300_dp, yg)
        syg = sqrt(yg)
        rinv = 1.0_dp / (sq * sh2 * syg)                 ! shared division
        inv_sq = sh2 * syg * rinv
        inv_sh2 = sq * syg * rinv
        inv_syg = sq * sh2 * rinv
        inv_z = inv_sq * inv_sq
        Cg = 2.0_dp * sh2 * sh2 * inv_z                  ! (1 - cos sqrt(g))/g
        sqC = sqrt2 * sh2 * inv_sq                       ! sqrt(C)
        invC = 0.5_dp * zg * inv_sh2 * inv_sh2         ! 1/C
        ! series below z = 0.5 avoids the sqrt(g) - sin(sqrt(g)) cancellation (both computed: no branch)
        Sg = merge(1.0_dp/6.0_dp + zg*(-1.0_dp/120.0_dp + zg*(1.0_dp/5040.0_dp &
                   + zg*(-1.0_dp/362880.0_dp + zg*(1.0_dp/39916800.0_dp &
                   + zg*(-1.0_dp/6227020800.0_dp + zg/1307674368000.0_dp))))), &
                   (sq - 2.0_dp * sh2 * ch2) * inv_z * inv_sq, zg < 0.5_dp)
        x3 = yg * syg * invC * (sq * inv_sh2 / sqrt2)   ! (y/C)^(3/2) = y sqrt(y) / (C sqrt(C))
        Phi = x3 * Sg + Ap * syg - mt
        dPhi = x3 * ((2.0_dp*Cg*Cg - 3.0_dp*Sg) * 0.25_dp * inv_z * invC + 0.75_dp * Sg * Sg * invC) &
             + 0.125_dp * Ap * (3.0_dp * Sg * syg * invC + Ap * sqC * inv_syg)
        dz = -Phi / dPhi
        dz = merge(dz, 0.0_dp, dz == dz .and. zg + dz > 0.0_dp .and. zg + dz < z_max)
        ! y and sqrt(y) at z* = g + dz to first order: dy/dz = sqrt(2) A sin(sqrt(g)/2) / (4 sqrt(g))
        dy = 0.25_dp * sqrt2 * Ap * sh2 * inv_sq * dz
        y = max(1.0e-300_dp, yg + dy)
        sy = syg + 0.5_dp * dy * inv_syg
    end subroutine polish

    subroutine solve_block(r1v, r2v, tof, sqrt_mu, inv_sqrt_mu, H, v1v, v2v, nb)
        integer, intent(in) :: nb
        real(dp), intent(in)  :: r1v(3,nb), r2v(3,nb), tof(nb), sqrt_mu, inv_sqrt_mu, H(3)
        real(dp), intent(out) :: v1v(3,nb), v2v(3,nb)

        integer :: j, p, k2, is(NBLK), kt(NBLK)
        logical :: bad(NBLK)
        real(dp), dimension(NBLK) :: inv_r1, inv_r2, r12, mu_t, A, s2absA, y_end, tp, I1_num, I0_den, z
        real(dp) :: r1_mag, r2_mag, c1, c2, c3, h_dot, dot12, d12sq, y0, ry0, z_new
        real(dp) :: wr, wi, qr, qi, tt, ar, ai, br, bi, cr, ci, dr, di, inv_d
        real(dp) :: xg, vg, fx, fy
        real(dp), dimension(NBLK) :: dz_p, y_p, sy_p
        integer :: ig, jg
        real(dp) :: sh, ch, sq, sh2, ch2, Cg, Sg, yg, syg, x3, Phi, dPhi, xa
        real(dp) :: rinv, inv_sq, inv_sh2, inv_syg, inv_z, sqC, invC, dy
        real(dp) :: ssq, csq, S_ser, S_cf, dzv, znew, ynew
        real(dp) :: f, gt, dg, inv_g

        ! Geometry
        do p = 1, nb
            r1_mag = sqrt(r1v(1,p)*r1v(1,p) + r1v(2,p)*r1v(2,p) + r1v(3,p)*r1v(3,p))
            r2_mag = sqrt(r2v(1,p)*r2v(1,p) + r2v(2,p)*r2v(2,p) + r2v(3,p)*r2v(3,p))
            inv_r1(p) = 1.0_dp / r1_mag
            inv_r2(p) = 1.0_dp / r2_mag
            r12(p) = r1_mag + r2_mag
            mu_t(p) = tof(p) * sqrt_mu

            c1 = r1v(2,p)*r2v(3,p) - r1v(3,p)*r2v(2,p)
            c2 = r1v(3,p)*r2v(1,p) - r1v(1,p)*r2v(3,p)
            c3 = r1v(1,p)*r2v(2,p) - r1v(2,p)*r2v(1,p)
            h_dot = c1*H(1) + c2*H(2) + c3*H(3)
            dot12 = r1v(1,p)*r2v(1,p) + r1v(2,p)*r2v(2,p) + r1v(3,p)*r2v(3,p)
            ! r1*r2 + r1.r2 cancels for transfer angles near 180 deg; there use the equal form
            ! |r1 x r2|^2 / (r1*r2 - r1.r2), whose denominator is ~2*r1*r2
            A(p) = r1_mag*r2_mag + dot12
            if (dot12 < 0.0_dp) A(p) = (c1*c1 + c2*c2 + c3*c3) / (r1_mag*r2_mag - dot12)
            A(p) = sqrt(max(0.0_dp, A(p)))
            if (h_dot < 0.0_dp) A(p) = -A(p)

            ! y(z) = y_end + sqrt(2)*|A|*(1 -/+ cos(sqrt(z)/2)) without cancellation:
            ! r12^2 - 2*(r1*r2 + r1.r2) = |r1 - r2|^2, so r12 - sqrt(2)*|A| = |r1 - r2|^2/(r12 + sqrt(2)*|A|)
            d12sq = (r1v(1,p)-r2v(1,p))**2 + (r1v(2,p)-r2v(2,p))**2 + (r1v(3,p)-r2v(3,p))**2
            s2absA(p) = sqrt2 * abs(A(p))
            y_end(p) = d12sq / (r12(p) + s2absA(p))
            is(p) = 2
            if (A(p) > 0.0_dp) is(p) = 1
            y0 = r12(p) + s2absA(p)                             ! y(0) = r12 - sqrt(2)*A
            if (A(p) > 0.0_dp) y0 = y_end(p)
            ry0 = sqrt(y0)
            tp(p) = (c_tp * y0 * ry0 + A(p) * ry0) * inv_sqrt_mu
            bad(p) = tof(p) < tp(p) .or. abs(A(p)) < 1.0e-12_dp .or. y_end(p) < sqrt2 * tau_gap * r12(p)
        end do

        ! Guess: z from the guess table (bilinear); near z = 0 the Taylor step z = (dt - tp)/T'(0),
        ! sqrt(mu) T'(0) = A sqrt(y0)/8 + sqrt(2) y0^(3/2)/40 + sqrt(2) A^2/(16 sqrt(y0)), y0 = y(0)
        do p = 1, nb
            xg = sqrt(y_end(p) / r12(p)) * (Nx_g - 1)
            vg = sqrt(max(0.0_dp, min(1.0_dp, tp(p) / tof(p)))) * (Nv_g - 1)
            ig = max(0, min(Nx_g - 2, int(xg)))
            jg = max(0, min(Nv_g - 2, int(vg)))
            fx = xg - ig
            fy = vg - jg
            z(p) = (1.0_dp - fy) * ((1.0_dp - fx) * zg_tab(ig, jg, is(p)) + fx * zg_tab(ig + 1, jg, is(p))) &
                 + fy * ((1.0_dp - fx) * zg_tab(ig, jg + 1, is(p)) + fx * zg_tab(ig + 1, jg + 1, is(p)))
            if (z(p) < z_thr) then
                y0 = r12(p) + s2absA(p)
                if (A(p) > 0.0_dp) y0 = y_end(p)
                ry0 = sqrt(y0)
                z(p) = (mu_t(p) - sqrt_mu * tp(p)) &
                     / (0.125_dp * A(p) * ry0 + 0.025_dp * sqrt2 * y0 * ry0 + 0.0625_dp * sqrt2 * A(p) * A(p) / ry0)
            end if
            if (.not. (z(p) > 0.0_dp)) z(p) = 1.0e-300_dp
            z(p) = min(z_max * (1.0_dp - 1.0e-16_dp), z(p))
        end do

        ! Contour: pick the small rectangle
        do p = 1, nb
            if (z(p) < z_thr) then
                kt(p) = id_end0 + min(j_hi, max(j_lo, exponent(z(p) / 1.5_dp * sqrt2) - 1)) - j_lo
            else if (z_max - z(p) < z_thr) then
                kt(p) = id_end4 + min(j_hi, max(j_lo, exponent((z_max - z(p)) / 1.5_dp * sqrt2) - 1)) - j_lo
            else
                kt(p) = min(N_tile, max(0, nint(z(p) * N_tile / z_max)))
            end if
        end do

        ! Contour: small rectangle, upper half only (num, den reuse I1_num, I0_den)
        I1_num = 0.0_dp
        I0_den = 0.0_dp
        do j = 1, NP
            do p = 1, nb
                wr = y_end(p) + s2absA(p) * tu_re(j, kt(p), is(p))
                wi = s2absA(p) * tu_im(j, kt(p), is(p))
                qr = sqrt(wr*wr + wi*wi) + wr                ! branch-free so the loop vectorizes
                tt = sqrt(2.0_dp * qr)
                qi = wi
                ar = wr*qr - wi*qi
                ai = wr*qi + wi*qr
                br = ar*tK_re(j, kt(p)) - ai*tK_im(j, kt(p))
                bi = ar*tK_im(j, kt(p)) + ai*tK_re(j, kt(p))
                cr = A(p)*qr - tt*mu_t(p)
                ci = A(p)*qi
                dr = br + cr*tq_re(j, kt(p)) - ci*tq_im(j, kt(p))
                di = bi + cr*tq_im(j, kt(p)) + ci*tq_re(j, kt(p))
                inv_d = tt / (dr*dr + di*di + 1.0e-300_dp)
                I1_num(p) = I1_num(p) + (tWz_im(j, kt(p))*dr - tWz_re(j, kt(p))*di) * inv_d
                I0_den(p) = I0_den(p) + (tW_im(j, kt(p))*dr - tW_re(j, kt(p))*di) * inv_d
            end do
        end do

        do p = 1, nb
            z_new = tile_c(kt(p)) + I1_num(p) / I0_den(p)
            if (z_new == z_new .and. z_new > tile_a(kt(p)) .and. z_new < tile_b(kt(p))) then
                z(p) = z_new
            else if (.not. bad(p) .and. z_new > 0.0_dp .and. z_new < z_max) then
                k2 = pick_tile(z_new)                       ! guess chose the wrong tile (rare): one retry
                if (k2 /= kt(p)) then
                    call top_sums(NP, tu_re(:, k2, is(p)), tu_im(:, k2, is(p)), tK_re(:, k2), tK_im(:, k2), &
                                  tq_re(:, k2), tq_im(:, k2), tW_re(:, k2), tW_im(:, k2), tWz_re(:, k2), tWz_im(:, k2), &
                                  y_end(p), s2absA(p), A(p), mu_t(p), I1_num(p), I0_den(p))
                    z_new = tile_c(k2) + I1_num(p) / I0_den(p)
                    if (z_new == z_new .and. z_new > tile_a(k2) .and. z_new < tile_b(k2)) then
                        z(p) = z_new
                        kt(p) = k2
                    end if
                end if
            end if

        end do

        ! Single polish step for the whole block: written out inline (no call, no table) so it vectorizes
        do p = 1, nb
            sq = sqrt(z(p))
            ! sin, cos of x = sqrt(g)/4 in [0, pi/2] by odd polynomials (no table, so the loop vectorizes);
            ! cos x = sin(pi/2 - x) with pi/2 in two parts keeps full relative accuracy near z = 4*pi^2
            xa = 0.25_dp * sq
            sh = sin_poly(xa)
            ch = sin_poly((pio2_hi - xa) + pio2_lo)
            sh2 = 2.0_dp * sh * ch                           ! sin(sqrt(g)/2)
            ch2 = (ch - sh) * (ch + sh)                      ! cos(sqrt(g)/2)
            ! both candidates are computed first, so each merge only selects finished values and the
            ! compiler can turn it into a vector blend (it will not speculate arithmetic inside a merge)
            ssq = sh * sh
            csq = ch * ch
            yg = y_end(p) + 2.0_dp * s2absA(p) * merge(ssq, csq, A(p) > 0.0_dp)   ! y(g) without cancellation
            yg = merge(yg, 1.0e-300_dp, yg > 1.0e-300_dp)      ! merge, not max: gfortran's max adds NaN jumps
            syg = sqrt(yg)
            rinv = 1.0_dp / (sq * sh2 * syg)                 ! shared division
            inv_sq = sh2 * syg * rinv
            inv_sh2 = sq * syg * rinv
            inv_syg = sq * sh2 * rinv
            inv_z = inv_sq * inv_sq
            Cg = 2.0_dp * sh2 * sh2 * inv_z                  ! (1 - cos sqrt(g))/g
            sqC = sqrt2 * sh2 * inv_sq                       ! sqrt(C)
            invC = 0.5_dp * z(p) * inv_sh2 * inv_sh2         ! 1/C
            ! series below z = 0.5 avoids the sqrt(g) - sin(sqrt(g)) cancellation (both computed: no branch)
            S_ser = 1.0_dp/6.0_dp + z(p)*(-1.0_dp/120.0_dp + z(p)*(1.0_dp/5040.0_dp &
                  + z(p)*(-1.0_dp/362880.0_dp + z(p)*(1.0_dp/39916800.0_dp &
                  + z(p)*(-1.0_dp/6227020800.0_dp + z(p)/1307674368000.0_dp)))))
            S_cf = (sq - 2.0_dp * sh2 * ch2) * inv_z * inv_sq
            Sg = merge(S_ser, S_cf, z(p) < 0.5_dp)
            x3 = yg * syg * invC * (sq * inv_sh2 / sqrt2)   ! (y/C)^(3/2) = y sqrt(y) / (C sqrt(C))
            Phi = x3 * Sg + A(p) * syg - mu_t(p)
            dPhi = x3 * ((2.0_dp*Cg*Cg - 3.0_dp*Sg) * 0.25_dp * inv_z * invC + 0.75_dp * Sg * Sg * invC) &
                 + 0.125_dp * A(p) * (3.0_dp * Sg * syg * invC + A(p) * sqC * inv_syg)
            dzv = -Phi / dPhi
            znew = z(p) + dzv
            dzv = merge(dzv, 0.0_dp, dzv == dzv)
            dzv = merge(dzv, 0.0_dp, znew > 0.0_dp)
            dzv = merge(dzv, 0.0_dp, znew < z_max)
            dz_p(p) = dzv
            ! y and sqrt(y) at z* = g + dz to first order: dy/dz = sqrt(2) A sin(sqrt(g)/2) / (4 sqrt(g))
            dy = 0.25_dp * sqrt2 * A(p) * sh2 * inv_sq * dzv
            ynew = yg + dy
            y_p(p) = merge(ynew, 1.0e-300_dp, ynew > 1.0e-300_dp)
            sy_p(p) = syg + 0.5_dp * dy * inv_syg
        end do

        do p = 1, nb
            ! A polish step far larger than the contour error means the tile held no root and its
            ! estimate fell inside it by chance (rare): solve on the tile containing z + dz, polish there
            if (.not. bad(p) .and. abs(dz_p(p)) > 1.0e-3_dp * (tile_b(kt(p)) - tile_a(kt(p)))) then
                k2 = pick_tile(z(p) + dz_p(p))
                if (k2 /= kt(p)) then
                    call top_sums(NP, tu_re(:, k2, is(p)), tu_im(:, k2, is(p)), tK_re(:, k2), tK_im(:, k2), &
                                  tq_re(:, k2), tq_im(:, k2), tW_re(:, k2), tW_im(:, k2), tWz_re(:, k2), tWz_im(:, k2), &
                                  y_end(p), s2absA(p), A(p), mu_t(p), I1_num(p), I0_den(p))
                    z_new = tile_c(k2) + I1_num(p) / I0_den(p)
                    if (z_new == z_new .and. z_new > tile_a(k2) .and. z_new < tile_b(k2)) then
                        z(p) = z_new
                        call polish(z(p), y_end(p), s2absA(p), A(p), mu_t(p), dz_p(p), y_p(p), sy_p(p))
                    end if
                end if
            end if

        end do

        ! Terminal velocities
        do p = 1, nb
            f = 1.0_dp - y_p(p) * inv_r1(p)
            gt = A(p) * sy_p(p) * inv_sqrt_mu
            dg = 1.0_dp - y_p(p) * inv_r2(p)
            inv_g = 1.0_dp / gt
            inv_g = merge(nan_val, inv_g, bad(p))
            inv_g = merge(nan_val, inv_g, abs(gt) < 1.0e-12_dp)
            v1v(1,p) = (r2v(1,p) - f*r1v(1,p)) * inv_g       ! components written out: no inner loop
            v1v(2,p) = (r2v(2,p) - f*r1v(2,p)) * inv_g
            v1v(3,p) = (r2v(3,p) - f*r1v(3,p)) * inv_g
            v2v(1,p) = (dg*r2v(1,p) - r1v(1,p)) * inv_g
            v2v(2,p) = (dg*r2v(2,p) - r1v(2,p)) * inv_g
            v2v(3,p) = (dg*r2v(3,p) - r1v(3,p)) * inv_g
        end do
    end subroutine solve_block
end module rect_mod_fast
