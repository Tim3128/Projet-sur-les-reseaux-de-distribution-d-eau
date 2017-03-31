exec('Wolfe_Skel.sci')

function [fopt,xopt,gopt]=Polack_Ribiere(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de Polack_Ribiere à pas variables                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de la méthode de Polack_Ribière a pas variables";
   labels = ["Nombre maximal d''iterations";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1);
   default = ["5000";"0.000001"];
   [ok,iter,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   xk_1 = x;
   xk = x;
   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [Fk,Gk] = Oracle(xk,ind);

//    - test de convergence

      if norm(Gk) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente
      if k = 1 then
          Dk = -Gk;
      else
          Betak = Gk'*(Gk-Gk_1)/(norm(Gk_1)^2)
          Dk = -Gk + Betak * Dk_1;
      end


//    - calcul de la longueur du pas de gradient
    
      //deltak = (F+4)/(k^(5/3));
      //alphak0 = (-2)*deltak/(G'*D);
      [alpha,ok] = Wolfe(1,xk,Dk,Oracle);

//    - mise a jour des variables
      Gk_1 = Gk;
      Dk_1 = Dk;
      xk_1 = xk;
      xk = xk + (alpha*Dk);

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(Gk)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; Fk ];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = Fk;
   xopt = xk;
   gopt = Gk;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de gradient a pas variables')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction
