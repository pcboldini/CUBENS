  rhoa  = r(i,j,k)
  ua    = u(i,j,k)
  va    = v(i,j,k)
  wa    = w(i,j,k)
  pa    = p(i,j,k)
  ea    = e(i,j,k)
  za    = zp(k)
  xa    = xp(i)

  !$acc loop seq
  do s=1,nStencilConv
    rhoip = rhoa+r(i+s,j,k)
    rhoim = rhoa+r(i-s,j,k)
    rhojp = rhoa+r(i,j+s,k)
    rhojm = rhoa+r(i,j-s,k)
    rhokp = rhoa+r(i,j,k+s)
    rhokm = rhoa+r(i,j,k-s)
    uip = ua+u(i+s,j,k)
    uim = ua+u(i-s,j,k)
    ujp = ua+u(i,j+s,k)
    ujm = ua+u(i,j-s,k)
    ukp = ua+u(i,j,k+s)
    ukm = ua+u(i,j,k-s)
    vip = va+v(i+s,j,k)
    vim = va+v(i-s,j,k)
    vjp = va+v(i,j+s,k)
    vjm = va+v(i,j-s,k)
    vkp = va+v(i,j,k+s)
    vkm = va+v(i,j,k-s)
    wip = wa+w(i+s,j,k)
    wim = wa+w(i-s,j,k)
    wjp = wa+w(i,j+s,k)
    wjm = wa+w(i,j-s,k)
    wkp = wa+w(i,j,k+s)
    wkm = wa+w(i,j,k-s)
    eip = ea+e(i+s,j,k)
    eim = ea+e(i-s,j,k)
    ejp = ea+e(i,j+s,k)
    ejm = ea+e(i,j-s,k)
    ekp = ea+e(i,j,k+s)
    ekm = ea+e(i,j,k-s)
    rhs_r(i,j,k) = rhs_r(i,j,k) - xa*0.5*conv_ddx(s)*( rhoip*uip - rhoim*uim ) &
                                -       0.5*conv_ddy(s)*( rhojp*vjp - rhojm*vjm ) &
                                - za*0.5*conv_ddz(s)*( rhokp*wkp - rhokm*wkm ) 
    rhs_u(i,j,k) = rhs_u(i,j,k) - xa*0.25*conv_ddx(s)*( rhoip*uip*uip - rhoim*uim*uim ) &
                                -       0.25*conv_ddy(s)*( rhojp*vjp*ujp - rhojm*vjm*ujm ) &
                                - za*0.25*conv_ddz(s)*( rhokp*wkp*ukp - rhokm*wkm*ukm ) 
    rhs_v(i,j,k) = rhs_v(i,j,k) - xa*0.25*conv_ddx(s)*( rhoip*uip*vip - rhoim*uim*vim ) &
                                -       0.25*conv_ddy(s)*( rhojp*vjp*vjp - rhojm*vjm*vjm ) &
                                - za*0.25*conv_ddz(s)*( rhokp*wkp*vkp - rhokm*wkm*vkm ) 
    rhs_w(i,j,k) = rhs_w(i,j,k) - xa*0.25*conv_ddx(s)*( rhoip*uip*wip - rhoim*uim*wim ) &
                                -       0.25*conv_ddy(s)*( rhojp*vjp*wjp - rhojm*vjm*wjm ) &
                                - za*0.25*conv_ddz(s)*( rhokp*wkp*wkp - rhokm*wkm*wkm ) 
    rhs_e(i,j,k) = rhs_e(i,j,k) - xa*0.25*conv_ddx(s)*(rhoip*uip*(eip)-rhoim*uim*(eim) &
                                             +rhoip*uip*(ua*u(i+s,j,k)+va*v(i+s,j,k)+wa*w(i+s,j,k)) &
                                             -rhoim*uim*(ua*u(i-s,j,k)+va*v(i-s,j,k)+wa*w(i-s,j,k))) &
                                - xa*conv_ddx(s)*( (ua*p(i+s,j,k) + pa*u(i+s,j,k)) &
                                            - (ua*p(i-s,j,k) + pa*u(i-s,j,k)) ) &
                                ! 
                                -      0.25*conv_ddy(s)*(rhojp*vjp*(ejp)-rhojm*vjm*(ejm) &
                                             +rhojp*vjp*(ua*u(i,j+s,k)+va*v(i,j+s,k)+wa*w(i,j+s,k)) &
                                             -rhojm*vjm*(ua*u(i,j-s,k)+va*v(i,j-s,k)+wa*w(i,j-s,k))) &
                                -      conv_ddy(s)*( (va*p(i,j+s,k) + pa*v(i,j+s,k)) &
                                            - (va*p(i,j-s,k) + pa*v(i,j-s,k)) ) &
                                ! 
                                - za*0.25*conv_ddz(s)*(rhokp*wkp*(ekp)-rhokm*wkm*(ekm) &
                                             +rhokp*wkp*(ua*u(i,j,k+s)+va*v(i,j,k+s)+wa*w(i,j,k+s)) &
                                             -rhokm*wkm*(ua*u(i,j,k-s)+va*v(i,j,k-s)+wa*w(i,j,k-s))) &
                                - za*conv_ddz(s)*( (wa*p(i,j,k+s) + pa*w(i,j,k+s)) &
                                            - (w(i,j,k)*p(i,j,k-s) + pa*w(i,j,k-s)) )
    rhs_u(i,j,k) = rhs_u(i,j,k) - xa*conv_ddx(s)*(p(i+s,j,k)-p(i-s,j,k)) 
    rhs_v(i,j,k) = rhs_v(i,j,k) -       conv_ddy(s)*(p(i,j+s,k)-p(i,j-s,k)) 
    rhs_w(i,j,k) = rhs_w(i,j,k) - za*conv_ddz(s)*(p(i,j,k+s)-p(i,j,k-s)) 
  enddo
