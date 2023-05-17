function llaPos = ecef2lla(ecefPos)

    [a, f] = fusion.internal.frames.wgs84ModelParams(class(ecefPos));

    x = ecefPos(:,1);
    y = ecefPos(:,2);
    z = ecefPos(:,3);

    rho = hypot(x, y);

    lon = atan2d(y, x);

    [lat, alt] = map.geodesy.internal.cylindrical2geodetic(rho, z, a, f, true);

    llaPos = [lat lon alt];
end