function glue_full_torsion(f,g,roots_f,roots_g)
     Qx := Parent(f);
     g := Qx ! g;
     x := Qx.1;    

    //roots_f := [ee[1] : ee in Roots(f)];
    //roots_g := [ee[1] : ee in Roots(g)];

    assert #(roots_f) eq 3;
    assert #(roots_g) eq 3;
    a1,a2,a3 := Explode(roots_f);
    b1,b2,b3 := Explode(roots_g);
    

    A1 := (a3-a2)^2/(b3-b2) + (a2-a1)^2/(b2-b1) + (a1-a3)^2/(b1-b3);
    B1 := (b3-b2)^2/(a3-a2) + (b2-b1)^2/(a2-a1) + (b1-b3)^2/(a1-a3);
    A2 := a1*(b3-b2) + a2*(b1-b3) + a3*(b2-b1);
    B2 := b1*(a3-a2) + b2*(a1-a3) + b3*(a2-a1);
    delta_f := (a1-a2)^2*(a2-a3)^2*(a1-a3)^2;
    delta_g := (b1-b2)^2*(b2-b3)^2*(b1-b3)^2;
    AA := delta_g*A1/A2;
    BB := delta_f*B1/B2;
    h := - (AA*(a2-a1)*(a1-a3)*x^2+BB*(b2-b1)*(b1-b3))*\
          (AA*(a3-a2)*(a2-a1)*x^2+BB*(b3-b2)*(b2-b1))*\
          (AA*(a1-a3)*(a3-a2)*x^2+BB*(b1-b3)*(b3-b2));
    return h;
end function;



