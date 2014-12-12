function h = drest(r)
    %h = -2*R;
   h0 = 1.985;
   R = 1;
   h = -r.*(((h0^2+R^2)/(2*h0))^2-r.^2).^(-1/2);
%     h = -2*r;
end