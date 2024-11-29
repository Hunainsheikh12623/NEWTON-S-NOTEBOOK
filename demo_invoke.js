function toggleMenu() {
    const navbar = document.querySelector('.navbar');
    navbar.classList.toggle('active');
}

// class 9th:

// chapter 01 Kinematics:

function calcfinalvelocity() {
    let u = parseFloat(document.querySelector("#vi_finalvelocity").value);
    let a = parseFloat(document.querySelector("#a_finalvelocity").value);
    let t = parseFloat(document.querySelector("#t_finalvelocity").value);

    if (isNaN(u) || isNaN(a) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const finalvelocity = Module.ccall(
        "finalvelocity",
        'number',
        ['number','number','number'],
        [u,a,t]
    )

       
    console.log(finalvelocity);
    document.querySelector("#finding_finalvelocity").innerHTML = "final velocity is " + finalvelocity + " m/s";

        alert("vf = u + at \nvf = " + u + " + " + a + " x " + t + "\nvf = " + finalvelocity + " m/s");

};

function calcdisplacement() {
    let u = parseFloat(document.querySelector("#vi_displacement").value);
    let a = parseFloat(document.querySelector("#a_displacement").value);
    let t = parseFloat(document.querySelector("#t_displacement").value);

    if (isNaN(u) || isNaN(a) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const displacement = Module.ccall(
        "displacement",
        'number',
        ['number','number','number'],
        [u,a,t]
    )

    console.log(displacement);
    document.querySelector("#finding_displacement").innerHTML = "displacement is " + displacement + " m";

    alert("s = ut + 0.5at^2 \ns = " + u + " x " + t + " + 0.5 x " + a + " x " + t + "^2 \ns = " + displacement + " m");
};

function calcacceleration() {
    let v = parseFloat(document.querySelector("#vf_acceleration").value);
    let u = parseFloat(document.querySelector("#vi_acceleration").value);
    let s = parseFloat(document.querySelector("#s_acceleration").value);

    if (isNaN(v) || isNaN(u) || isNaN(s)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const acceleration = Module.ccall(
        "acceleration",
        'number',
        ['number','number','number'],
        [v,u,s]
    )

    console.log(acceleration);
    document.querySelector("#finding_acceleration").innerHTML = "acceleration is " + acceleration + " m/s^2";

    alert("a = (v^2 - u^2) / (2s) \na = (" + v + "^2 - " + u + "^2) / (2 x " + s + ") \na = " + acceleration + " m/s^2");
};

// motion under gravity:

function calcfinalvelocity_gravity() {
    let u = parseFloat(document.querySelector("#vi_gravity").value);
    let t = parseFloat(document.querySelector("#t_gravity").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const finalvelocity = Module.ccall(
        "finalvelocity_gravity",
        'number',
        ['number','number'],
        [u,t]
    );

    console.log(finalvelocity);
    document.querySelector("#finding_finalvelocity_gravity").innerHTML = "final velocity is " + finalvelocity + " m/s";

    alert("v = u + gt \nv = " + u + " + " + 9.8 + " x " + t + "\nv = " + finalvelocity + " m/s");
    
};

function calcheight_gravity() {
    let v = parseFloat(document.querySelector("#vf_gravity").value);

    if (isNaN(v)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const height = Module.ccall(
        "displacement_gravity",
        'number',
        ['number'],
        [v]
    );

    console.log(height);
    document.querySelector("#finding_height_gravity").innerHTML = "height is " + height + " m";

    alert("h = ut + 0.5gt^2 \nh = 0 + 0.5 x 9.8 x " + v + "^2 \nh = " + height + " m");
};

function calctime_gravity() {
    let h = parseFloat(document.querySelector("#h_gravity").value);

    if (isNaN(h)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const time = Module.ccall(
        "time_gravity",
        'number',
        ['number'],
        [h]
    );

    console.log(time);
    document.querySelector("#finding_time_gravity").innerHTML = "time is " + time + " s";
    document.querySelector("#finding_total_time_gravity").innerHTML = "total time is " + time*2 + " s";

    alert("h = ut + 0.5gt^2 \n" + h + " = 0 + 0.5 x 9.8 x t^2 \nt = " + "sqrt" + "(" + 2*h/9.8 + ")"+ "\nt = " + time + " s\ntotal = " + time*2 + " s");
};


// chapter 2: Dynamics:


function calcforce() {
    let m = parseFloat(document.querySelector("#m_force").value);
    let a = parseFloat(document.querySelector("#a_force").value);

    if (isNaN(m) || isNaN(a)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const force = Module.ccall(
        "force",
        'number',
        ['number','number'],
        [m,a]
    )

    console.log(force);
    document.querySelector("#finding_force").innerHTML = "force is " + force + " N";

        alert("F = ma \n" + "F = " + m + " x " + a + " \nF = " + force +" N");

};

function calcmass() {
    let f = parseFloat(document.querySelector("#f_mass").value);
    let a = parseFloat(document.querySelector("#a_mass").value);

    if (isNaN(f) || isNaN(a)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const mass = Module.ccall(
        "mass",
        'number',
        ['number','number'],
        [f,a]
    )

    console.log(mass);
    document.querySelector("#finding_mass").innerHTML = "mass is " + mass + " kg";

        alert("m = F/a \n" + "m = " + f + " / " + a + " \nm = " + mass +" kg");
};

function calcacceleration_force() {
    let f = parseFloat(document.querySelector("#f_acceleration").value);
    let m = parseFloat(document.querySelector("#m_acceleration").value);

    if (isNaN(f) || isNaN(m)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const acceleration = Module.ccall(
        "acceleration_force",
        'number',
        ['number','number'],
        [f,m]
    )

    console.log(acceleration);
    document.querySelector("#finding_acceleration_force").innerHTML = "acceleration is " + acceleration + " m/s^2";

        alert("a = F/m \n" + "a = " + f + " / " + m + " \na = " + acceleration +" m/s^2");
};

function calcmomentum() {
    let m = parseFloat(document.querySelector("#m_momentum").value);
    let v = parseFloat(document.querySelector("#v_momentum").value);

    if (isNaN(m) || isNaN(v)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const momentum = Module.ccall(
        "momentum",
        'number',
        ['number','number'],
        [m,v]
    );

    console.log(momentum);
    document.querySelector("#finding_momentum").innerHTML = "momentum is " + momentum + " kg*m/s";

        alert("p = m x v \n" + "p = " + m + " x " + v + " \np = " + momentum +" kg*m/s");
};


// chapter 3: Work, power and energy :

function calcwork() {
    let f = parseFloat(document.querySelector("#f_work").value);
    let d = parseFloat(document.querySelector("#d_work").value);

    if (isNaN(f) || isNaN(d)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const work = Module.ccall(
        "work",
        'number',
        ['number','number'],
        [f,d]
    )

    console.log(work);
    document.querySelector("#finding_work").innerHTML = "work is " + work + " Joules";

        alert("W = F x d \n" + "W = " + f + " x " + d + " \nW = " + work +" Joules");
};

function calcpower() {
    let w = parseFloat(document.querySelector("#w_power").value);
    let t = parseFloat(document.querySelector("#t_power").value);

    if (isNaN(w) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const power = Module.ccall(
        "power",
        'number',
        ['number','number'],
        [w,t]
    )

    console.log(power);
    document.querySelector("#finding_power").innerHTML = "power is " + power + " Watt";

        alert("P = W / t \n" + "P = " + w + " / " + t + " \nP = " + power +" Watt");
};

function calcenergy() {
    let m = parseFloat(document.querySelector("#m_energy").value);
    let v = parseFloat(document.querySelector("#v_energy").value);

    if (isNaN(m) || isNaN(v)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const energy = Module.ccall(
        "energy",
        'number',
        ['number','number'],
        [m,v]
    );

    console.log(energy);
    document.querySelector("#finding_energy").innerHTML = "energy is " + energy + " Joules";

        alert("E = 1/2 x m x v^2 \n"+ "K.E = 1/2 x "+ m + " x " + v + "^2 \nK.E = " + energy + " Joules")
};

function calc_p_energy() {
    let m = parseFloat(document.querySelector("#m_p_energy").value);
    let h = parseFloat(document.querySelector("#h_p_energy").value);

    if (isNaN(m) || isNaN(h)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const p_energy = Module.ccall(
        "p_energy",
        'number',
        ['number','number'],
        [m,h]
    );

    console.log(p_energy);
    document.querySelector("#finding_p_energy").innerHTML = "potential energy is " + p_energy + " Joules";

        alert("P.E = mgh \n"+ "P.E = "+ m + " x 9.8 x "  + h + "\nU = " + p_energy + " Joules")

};

// chapter 04 Gravitation:

function calc_force_gravitation() {
    let m1 = parseFloat(document.querySelector("#m1_gravitation").value);
    let e1 = parseFloat(document.querySelector("#m1_exp_gravitation").value);
    let m2 = parseFloat(document.querySelector("#m2_gravitation").value);
    let e2 = parseFloat(document.querySelector("#m2_exp_gravitation").value);
    let d  = parseFloat(document.querySelector("#d_gravitation").value);
    let e3 = parseFloat(document.querySelector("#d_exp_gravitation").value);

    if (d == 0) {
        alert("Distance cannot be zero.");
        return;
    }

    if (isNaN(m1) || isNaN(m2) || isNaN(d) || isNaN(e1) || isNaN(e2) || isNaN(e3)) {
        alert("Please enter valid numerical values.");
    }
        const force_gravitation = Module.ccall(
            "force_gravitation",
            'number',
            ['number','number','number','number','number','number'],
            [m1,e1,m2,e2,d,e3]
        );


        console.log(force_gravitation);
        document.querySelector("#finding_force_gravitation").innerHTML = "force of gravitation is " + force_gravitation + " N";


        alert("F = G x m1 x m2 / d^2 \n"+
              "F = 6.674 x 10^-11 x "+ m1 + "x" + "10^" + e1 + " x " + m2 + " x 10^" + e2 + " / " + "(" + d + " x " + "10^" + e3 + ")^2\n"+
              "F = " + "(" + "6.674" + "x" +  m1 + "X" + m2 + ")" + " / " + "(" + d + "^2" + ")" + "x 10^(" + "-11 + " + e1 + " + " + e2 + " - " + "2 x" + e3 + ")" + "\n"+
              "F = " + 6.674*m1*m2 + "/" + d*d + " x 10^" + (-11 + e1 + e2 - 2*e3) + "\n"+
              "F = " + force_gravitation + " N"
            );
};

function calc_m1_gravitation() {
    let m2 = parseFloat(document.querySelector("#m2_m1_gravitation").value);
    let e1 = parseFloat(document.querySelector("#m2_exp_m1_gravitation").value);
    let d  = parseFloat(document.querySelector("#d_m1_gravitation").value);
    let e2 = parseFloat(document.querySelector("#d_exp_m1_gravitation").value);
    let f  = parseFloat(document.querySelector("#f_m1_gravitation").value);
    let e3 = parseFloat(document.querySelector("#f_exp_m1_gravitation").value);

    if (m2 == 0) {
        alert("m2 cannot be zero");
    };

    if (isNaN(m2) || isNaN(e1) || isNaN(d) || isNaN(e2) || isNaN(f) || isNaN(e3)) {
        alert("Please enter valid numerical values.");
    };

    const m1 = Module.ccall(
        'm1_gravitation',
        'number',
        ['number','number','number','number','number','number'],
        [m2,e1,d,e2,f,e3]
    );

    console.log(m1);
    document.querySelector("#finding_m1_gravitation").innerHTML = "m1 is " + m1 + " kg";
};

function calc_d_gravitation() {
    let m1 = parseFloat(document.querySelector("#m1_d_gravitation").value);
    let e1 = parseFloat(document.querySelector("#m1_exp_d_gravitation").value);
    let m2 = parseFloat(document.querySelector("#m2_d_gravitation").value);
    let e2 = parseFloat(document.querySelector("#m2_exp_d_gravitation").value);
    let f  = parseFloat(document.querySelector("#f_d_gravitation").value);
    let e3 = parseFloat(document.querySelector("#f_exp_d_gravitation").value);

    if (f == 0) {
        alert("f cannot be zero");
    };

    if (isNaN(m1) || isNaN(e1) || isNaN(m2) || isNaN(e2) || isNaN(f) || isNaN(e3)) {
        alert("Please enter valid numerical values.");
    };

    const d = Module.ccall(
        'distance_gravitation',
        'number',
        ['number','number','number','number','number','number'],
        [m1,e1,m2,e2,f,e3]
    );

    console.log(d);
    document.querySelector("#finding_d_gravitation").innerHTML = "d is " + d + " m";
};




// class 10th:

// chapter 01: Electricity:

function calcvoltage() {
    let I = parseFloat(document.querySelector("#I_voltage").value);
    let R = parseFloat(document.querySelector("#R_voltage").value);

    if (isNaN(I) || isNaN(R)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const voltage = Module.ccall(
        "voltage",
        'number',
        ['number','number'],
        [I,R]
    );

    alert("V = IR\n"+"V ="+ I + " x " + R + "\n" + "V = " + voltage + " volts")

    console.log(voltage);
    document.querySelector("#finding_voltage").innerHTML = "voltage is " + voltage + " volts";
};

function calccurrent() {
    let V = parseFloat(document.querySelector("#V_current").value);
    let R = parseFloat(document.querySelector("#R_current").value);

    if (isNaN(V) || isNaN(R)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const current = Module.ccall(
        "current",
        'number',
        ['number','number'],
        [V,R]
    );

    alert("I = V/R\n"+"I ="+ V + " / " + R + "\n" + "I = " + current + " amperes")

    console.log(current);
    document.querySelector("#finding_current").innerHTML = "current is " + current + " amperes";

};

function calcresistance() {
    let V = parseFloat(document.querySelector("#V_resistance").value);
    let I = parseFloat(document.querySelector("#I_resistance").value);

    if (isNaN(V) || isNaN(I)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const resistance = Module.ccall(
        "resistance",
        'number',
        ['number','number'],
        [V,I]
    );

    alert("R = V/I\n"+"R ="+ V + " / " + I + "\n" + "R = " + resistance + " ohms")
    
    console.log(resistance);
    document.querySelector("#finding_resistance").innerHTML = "resistance is " + resistance + " ohms";

};

function calcEpower() {
    let V = parseFloat(document.querySelector("#V_Epower").value);
    let I = parseFloat(document.querySelector("#I_Epower").value);
    let R = parseFloat(document.querySelector("#R_Epower").value);

    if (isNaN(V) || isNaN(I) || isNaN(R)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const Epower = Module.ccall(
        "electricalpower",
        'number',
        ['number','number', 'number'],
        [V, I, R]
    );

    if(R === 0) {
        alert ("P = VI\n"+"P = "+ V +" x "+ I +"\n"+"P = "+ V*I +" watts");
    }
    else if( I === 0) {
        alert ("P = V^2/R\n"+"P = "+ V +"^2 / "+ R +"\n"+"P = "+ V*V/R +" watts");
    }
    else if( V === 0) {
        alert ("P = I^2*R\n"+"P = "+ I +"^2 x "+ R +"\n"+"P = "+ I*I*R +" watts");
    }
    else {
        alert ("Invalid Input");
    }

    console.log(Epower);
    document.querySelector("#finding_E_power").innerHTML = "Electrical power is " + Epower + " watts";
}

function calcEenergy() {
    let V = parseFloat(document.querySelector("#V_Eenergy").value);
    let I = parseFloat(document.querySelector("#I_Eenergy").value);
    let R = parseFloat(document.querySelector("#R_Eenergy").value);
    let T = parseFloat(document.querySelector("#t_Eenergy").value);

    if (isNaN(V) || isNaN(I) || isNaN(R) || isNaN(T)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const Eenergy = Module.ccall(
        "electricalenergy",
        'number',
        ['number','number', 'number', 'number'],
        [V, I, R, T]
    )

    if( R === 0) {
        alert ("E = VIt\n"+"E = "+ V +" x "+ I +" x "+ T +"\n"+"E = "+ V*I*T +" joules");
    }
    else if( I === 0) {
        alert ("E = V^2 * t /R\n"+"E = "+ V +"^2 x"+ T + " / "+ R  +"\n"+"E = "+ V*V*T/R+" joules");
    }
    else if( V === 0) {
        alert ("E = I^2 * R * t\n"+"E = "+ I +"^2 x "+ R +" x "+ T +"\n"+"E = "+ I*I*R*T +" joules");
    }
    else{
        alert("Time will not be zero!")
    }

    console.log(Eenergy);
    document.querySelector("#finding_E_energy").innerHTML = "Electrical Energy = " + Eenergy + " joules";
}


// class 11th
// chapter 01: Motion In two Dimensions

function calc_xcompo() {
    let u = parseFloat(document.querySelector("#V_components").value);
    let t = (document.querySelector("#theta_components").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const xcompo = Module.ccall(
        'hori_velocity_comp',
        'number',
        ['number', 'number'],
        [u, theta_rad]
    )

    alert("Vx = Vcos(theta)\nVx = "+ u +"cos("+ t +")\nVx = "+ xcompo);
    console.log(xcompo);
    document.querySelector("#finding_xcompo").innerHTML = "Horizontal Velocity = " + xcompo + " m/sec";
}

function calc_ycompo() {
    let u = parseFloat(document.querySelector("#V_components").value);
    let t = (document.querySelector("#theta_components").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const ycompo = Module.ccall(
        'vert_velocity_comp',
        'number',
        ['number', 'number'],
        [u, theta_rad]
    )



    alert("Vy = Vsin(theta)\nVy = "+ u +"sin("+ t +")\nVy = "+ ycompo);
    console.log(ycompo);
    document.querySelector("#finding_ycompo").innerHTML = "Vertical Velocity = " + ycompo + " m/sec";
}


function calc_resultant_velocity() {
    let Vx = parseFloat(document.querySelector("#Vx_Resultant").value);
    let Vy = parseFloat(document.querySelector("#Vy_Resultant").value);

    if (isNaN(Vx) || isNaN(Vy)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Calculate resultant velocity using Pythagorean theorem
    const resultant_velocity = Module.ccall(
        'resultant_velocity',
        'number',
        ['number', 'number'],
        [Vx, Vy]
    )

    // Display the correct formula and result in the alert box
    alert("Resultant Velocity = sqrt(Vx^2 + Vy^2)\n" +
        "Resultant Velocity = sqrt(" + Vx + "^2 + " + Vy + "^2)\n" +
        "Resultant Velocity = sqrt(" + Math.pow(Vx, 2) + " + " + Math.pow(Vy, 2) + ")\n" +
        "Resultant Velocity = " + resultant_velocity);

    console.log(resultant_velocity);

    // Display the calculated result on the webpage
    document.querySelector("#finding_resultant_velocity").innerHTML = "Resultant Velocity = " + resultant_velocity + " m/s";
}

function calc_Horizontal_range() {
    let u = parseFloat(document.querySelector("#V_range").value);
    let t = (document.querySelector("#theta_range").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const Range = Module.ccall(
        'horizontal_range',
        'number',
        ['number', 'number'],
        [u, theta_rad]
    )



    alert("R = V² x sin2(theta) / g \n" +" R ="+ u + "^2" +" x sin("+ 2 * t +")/ 9.8 \n  = " + "R = " + Range);
    console.log(Range);
    document.querySelector("#finding_horizontal_range").innerHTML = "Horizontal Range = " + Range + " m";


};

function calc_Maximum_height() {
    let u = parseFloat(document.querySelector("#V_height").value);
    let t = (document.querySelector("#theta_height").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const height = Module.ccall(
        'maximum_height',
        'number',
        ['number', 'number'],
        [u, theta_rad]
    )



    alert("R = V² x sin(theta)² / 2g \n" +" h ="+ u + "^2" +" x sin(" + t + ")" + "^ 2" + "/ "+ 2 * 9.8 + "\n " + "h = " + height);
    console.log(height);
    document.querySelector("#finding_maximum_height").innerHTML = "Maximum Height = " + height + " m";


};

function calc_Total_time() {
    let u = parseFloat(document.querySelector("#V_flight").value);
    let t = (document.querySelector("#theta_flight").value);

    if (isNaN(u) || isNaN(t)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const time = Module.ccall(
        'time_of_flight',
        'number',
        ['number', 'number'],
        [u, theta_rad]
    )



    alert("T = 2V x sin(theta) / g \n" +" T ="+ " 2 x " + u +" x sin(" + t + ")" + "/ "+ 9.8 + "\n " + "T = " + time);
    console.log(time);
    document.querySelector("#finding_total_time").innerHTML = "Total time = " + time + " sec";


};

function calc_positionX() {
    let u = parseFloat(document.querySelector("#V_Xposition").value);
    let t = (document.querySelector("#theta_cos").value);
    let S = parseFloat(document.querySelector("#time_horizontal").value);

    if (isNaN(u) || isNaN(t) || isNaN(S)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const Xposition = Module.ccall(
        'position_x_t',
        'number',
        ['number', 'number', 'number'],
        [u, theta_rad, S]
    )



    alert("Px = V x cos(theta) x t  \n" +" Px =" + u +" x cos(" + t + ")" + " x "+ S + "\n " + "Px = " + Xposition);
    console.log(Xposition);
    document.querySelector("#finding_Xposition").innerHTML = "Position on X axis = " + Xposition + " m";


};

function calc_positionY() {
    let u = parseFloat(document.querySelector("#V_Yposition").value);
    let t = (document.querySelector("#theta_sin").value);
    let S = parseFloat(document.querySelector("#time_verticle").value);

    if (isNaN(u) || isNaN(t) || isNaN(S)) {
        alert("Please enter valid numerical values.");
        return;
    }

    // Convert angle from degrees to radians
    let theta_rad = t * Math.PI / 180;

    const Yposition = Module.ccall(
        'position_y_t',
        'number',
        ['number', 'number', 'number'],
        [u, theta_rad, S]
    )



    alert("Py = (V x sin(theta) x t) - (0.5 x g x t^2)) \n" +" Py =" + u +" x sin(" + t + ")" + " x "+ S + "-"+ " (0.5 x 9.8 x "+ S + "^2)" + "\n " + "Py = " + Yposition);
    console.log(Yposition);
    document.querySelector("#finding_Yposition").innerHTML = "Position on Y axis = " + Yposition + " m";


};

// chapter 2, 3, 4 : sound, oscillations, waves

function calc_wave_speed() {
    let f = parseFloat(document.querySelector("#wave_frequency").value);
    let lamde = parseFloat(document.querySelector("#wave_wavelength").value);

    if (isNaN(f) || isNaN(lamde)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const wave_speed = Module.ccall(
        'wave_speed',
        'number',
        ['number', 'number'],
        [f, lamde]
    );


    alert("V = f x lamda \n" + "Wave speed = " + f + " x " + lamde + "\n " + "Wave speed = " + wave_speed + " m/s");
    console.log(wave_speed);
    document.querySelector("#finding_wave_speed").innerHTML = "Wave speed = " + wave_speed + " m/s";
};

function calc_frequency_from_period() {
    let T = parseFloat(document.querySelector("#time_period").value);

    if (isNaN(T)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const frequency = Module.ccall(
        'frequency_from_period',
        'number',
        ['number'],
        [T]
    );

    alert("f = 1/T \n" + "f = " + 1 + " / " + T + "\n " + "f = " + frequency + " Hz");
    console.log(frequency);
    document.querySelector("#finding_frequency").innerHTML = "Frequency = " + frequency + " Hz";
};


function calc_kEshm() {
    let m = parseFloat(document.querySelector("#m_shm").value);
    let w = parseFloat(document.querySelector("#w_shm").value);
    let A = parseFloat(document.querySelector("#A_shm").value);
    let X = parseFloat(document.querySelector("#X_shm").value);

    if (isNaN(m) || isNaN(w) || isNaN(A) || isNaN(X)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const kineticEnergySHM = Module.ccall(
        'kinetic_energy_shm',
        'number',
        ['number', 'number', 'number', 'number'],
        [m, w, A, X]
    );

    alert(("KE = 1/2 m w^2 (A^2 - X^2)\n KE = 1/2 ") + m + " x " + w + "^2 ( " + A + "^2 - " + X + "^2 ) \n KE = " + kineticEnergySHM + " J")

    console.log(kineticEnergySHM);
    document.querySelector("#finding_kinetic_energy_shm").innerHTML = "Kinetic Energy = " + kineticEnergySHM + " J";
};

function calc_PEshm() {
    let m = parseFloat(document.querySelector("#mp_shm").value);
    let w = parseFloat(document.querySelector("#wp_shm").value);
    let X = parseFloat(document.querySelector("#Xp_shm").value);

    if (isNaN(m) || isNaN(w) || isNaN(X)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const potentialEnergySHM = Module.ccall(
        'potential_energy_shm',
        'number',
        ['number', 'number', 'number'],
        [m, w, X]
    );

    alert(("PE = 1/2 * m * ω² * x²\n PE= 1/2 ") + m + " x " + w + "^2 x " + X + "^2 ) \n PE = " + potentialEnergySHM + " J")

    console.log(potentialEnergySHM);
    document.querySelector("#finding_potential_energy_shm").innerHTML = "Kinetic Energy = " + potentialEnergySHM + " J";
};

function calc_wp() {
    let m = parseFloat(document.querySelector("#md_wp").value);
    let w = parseFloat(document.querySelector("#w_wp").value);
    let A = parseFloat(document.querySelector("#A_wp").value);
    let V = parseFloat(document.querySelector("#V_wp").value);


    if (isNaN(m) || isNaN(w) || isNaN(A) || isNaN(V)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const wavePower = Module.ccall(
        'wave_power',
        'number',
        ['number', 'number', 'number', 'number'],
        [m, w, A, V]
    );

    alert("WP = 1/2 * μ * ω² * A² * V\n WP= 1/2 " + m + " x " + w + "^2 x " + A + "^2 x " + V + "\n WP = " + wavePower + " J")

    console.log(wavePower);
    document.querySelector("#finding_wave_power").innerHTML = "Wave Power  = " + wavePower + " watt";
};

function calc_de_Tsource() {
    let fs = parseFloat(document.querySelector("#sf_tsource").value);
    let V = parseFloat(document.querySelector("#V_tsource").value);
    let Vo = parseFloat(document.querySelector("#Vo_tsource").value);

    if (isNaN(fs) || isNaN(V) || isNaN(Vo)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const de_Tsource = Module.ccall(
        'de_Tsource',
        'number',
        ['number', 'number', 'number'],
        [fs, V, Vo]
    )

    alert("f' = f (v + vo) / v \n f' = " + fs + " x (" + V + " + " + Vo + ") / " + V + "\n f' = " + de_Tsource + " Hz");
    console.log(de_Tsource);
    document.querySelector("#finding_de_tsource").innerHTML = "f' = " + de_Tsource + " Hz" 

};

function calc_de_Asource() {
    let fs = parseFloat(document.querySelector("#sf_asource").value);
    let V = parseFloat(document.querySelector("#V_asource").value);
    let Vo = parseFloat(document.querySelector("#Vo_asource").value);

    if (isNaN(fs) || isNaN(V) || isNaN(Vo)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const de_Asource = Module.ccall(
        'de_Asource',
        'number',
        ['number', 'number', 'number'],
        [fs, V, Vo]
    )

    alert("f' = f (v - vo) / v \n f' = " + fs + " x (" + V + " - " + Vo + ") / " + V + "\n f' = " + de_Asource + " Hz");
    console.log(de_Asource);
    document.querySelector("#finding_de_asource").innerHTML = "f' = " + de_Asource + " Hz" 

};

function calc_de_Tobserver() {
    let fs = parseFloat(document.querySelector("#sf_tobserver").value);
    let V = parseFloat(document.querySelector("#V_tobserver").value);
    let Vs = parseFloat(document.querySelector("#Vs_tobserver").value);

    if (isNaN(fs) || isNaN(V) || isNaN(Vs)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const de_Tobserver = Module.ccall(
        'de_Tobserver',
        'number',
        ['number', 'number', 'number'],
        [fs, V, Vs]
    )

    alert("f' = f V / (v - vs) \n f' = " + fs + " x " + V +  "/" + "(" + V + " - " + Vs +")  \n  f' = " + de_Tobserver + " Hz");
    console.log(de_Tobserver);
    document.querySelector("#finding_de_tobserver").innerHTML = "f' = " + de_Tobserver + " Hz" 

};

function calc_de_Aobserver() {
    let fs = parseFloat(document.querySelector("#sf_aobserver").value);
    let V = parseFloat(document.querySelector("#V_aobserver").value);
    let Vs = parseFloat(document.querySelector("#Vs_aobserver").value);

    if (isNaN(fs) || isNaN(V) || isNaN(Vs)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const de_Aobserver = Module.ccall(
        'de_Aobserver',
        'number',
        ['number', 'number', 'number'],
        [fs, V, Vs]
    )

    alert("f' = f V / (v + vs) \n f' = " + fs + " x " + V +  "/" + "(" + V + " + " + Vs +")  \n  f' = " + de_Aobserver + " Hz");
    console.log(de_Aobserver);
    document.querySelector("#finding_de_aobserver").innerHTML = "f' = " + de_Aobserver + " Hz" 

};

function calc_beat_frequency() {
    let f1 = parseFloat(document.querySelector("#f1_beat").value);
    let f2 = parseFloat(document.querySelector("#f2_beat").value);

    if (isNaN(f1) || isNaN(f2)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const beat_frequency = Module.ccall(
        'beat_frequency',
        'number',
        ['number', 'number'],
        [f1, f2]
    );

    alert("F = |f1 - f2| \n F = |"+ f1 + " - " + f2 + "| \n F = " + beat_frequency + " Hz");
    console.log(beat_frequency);
    document.querySelector("#finding_beat_frequency").innerHTML = "F = " + beat_frequency + " Hz"
};

function calc_standing_wave_wavelength() {
    let l = parseFloat(document.querySelector("#l_sw").value);
    let n = parseFloat(document.querySelector("#n_sw").value);

    if (isNaN(l) || isNaN(n)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const sw_wavelength = Module.ccall(
        'standing_wave_wavelength',
        'number',
        ['number', 'number'],
        [l, n]
    );

    alert("λ = 2l / n \n λ = 2 * " + l + " / " + n + "\n λ = " + sw_wavelength + " m");
    console.log(sw_wavelength);
    document.querySelector("#finding_standing_wave_wavelength").innerHTML = "λ = " + sw_wavelength + " m"
    console.log(sw_wavelength);
    document.querySelector("#finding_standing_wave_wavelength").innerHTML = "λ = " + sw_wavelength + " m"
};

function calc_resonance_frequency() {
    let K = parseFloat(document.querySelector("#kf").value);
    let m = parseFloat(document.querySelector("#mf").value);

    if (isNaN(K) || isNaN(m)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const res_frequency = Module.ccall(
        'resonance_frequency',
        'number',
        ['number', 'number'],
        [K, m]
    );

    alert("f = 1 / 2π √(k / m) \n"+ "f = 1 / 2π √(" + K + " / " + m + ") \n f = " + res_frequency + " Hz");
    console.log(res_frequency);
    document.querySelector("#finding_resonance_frequency").innerHTML = "f = " + res_frequency + " Hz"
};

function calc_sound_intensity() {
    let P = parseFloat(document.querySelector("#ps").value);
    let A = parseFloat(document.querySelector("#as").value);

    if (isNaN(P) || isNaN(A)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const sound_intensity = Module.ccall(
        'sound_intensity',
        'number',
        ['number', 'number'],
        [P, A]
    );

    alert("I = P / A \n I = "+ P + " / " + A + "\n I = " + sound_intensity + " W/m²");
    console.log(sound_intensity);
    document.querySelector("#finding_sound_intensity").innerHTML = "I = " + sound_intensity + " W/m²"
};

function calc_sound_level() {
    let Is = parseFloat(document.querySelector("#Is").value);
    let Io = parseFloat(document.querySelector("#Io").value);

    if (isNaN(Is) || isNaN(Io)) {
        alert("Please enter valid numerical values.");
        return;
    }

    const sound_level = Module.ccall(
        'sound_level',
        'number',
        ['number', 'number'],
        [Is, Io]
    );

    alert("L = 10 * log10(I / Io) \n L = 10 * log(" + Is + " / " + Io + ") \n L = " + sound_level + " dB");

    console.log(sound_level);
    document.querySelector("#finding_sound_level").innerHTML = "L = " + sound_level + " dB";
}

// Electric Field Calculation
function calc_electric_field() {
    const F = parseFloat(document.getElementById('force_electric').value);
    const q = parseFloat(document.getElementById('charge_electric').value);
    const r = parseFloat(document.getElementById('distance_field').value);

    if (isNaN(F) && isNaN(q) && isNaN(r)) {
        alert("Please enter valid values for Force, Charge, and Distance.");
    }

    const electricField = Module.ccall(
        'electric_field', 
        'number', 
        ['number', 'number', 'number'], 
        [F, q, r]
    );
    
        if(r === 0) {
            alert("E = F / q \n E = " + F + " / " + q + " \n E = " + electricField + " N/C" )
        }
        else if(F === 0) {
            alert("E = k x q / r^2 \n E = " + "9 x 10^9 x" + q + " / " + r + "^2 \n E = " + electricField + " N/C" )
        }

        console.log(electricField);
        document.querySelector("#finding_electric_field").innerHTML = "E = " + electricField + " N/C";
}

// Electric Potential Calculation
function calc_electric_potential() {
    const q = parseFloat(document.getElementById('charge_potential').value);
    const r = parseFloat(document.getElementById('distance_potential').value);

    if (isNaN(q) && isNaN(r)) {
        alert("Please enter valid values for Charge and Distance.");
    } 

    const electricPotential = Module.ccall(
        'electric_potential', 
        'number', 
        ['number', 'number'], 
        [q, r]
    );

    alert("V = K x q / r \n V = 9 x 10^9 x " + q + " / " + r + " \n V = " + electricPotential + " J/C")

    console.log(electricPotential);
    document.querySelector("#finding_electric_potential").innerHTML = "V = " + electricPotential + " J/C";
}

// Capacitance Calculation
function calc_capacitance() {
    const Q = parseFloat(document.getElementById('charge_capacitance').value);
    const V = parseFloat(document.getElementById('voltage_capacitance').value);
    
    const capacitance = Module.ccall(
        'capacitance', 
        'number', 
        ['number', 'number', 'number', 'number'], 
        [Q, V, 0, 0]
    );

    if (!isNaN(Q) && !isNaN(V)) {
       alert("C = Q / V \n C = " + Q + " / " + V + " \n C = " + capacitance + " F");
    }
    else {
        alert("Please enter valid values.");
    }

    console.log(capacitance);
    document.querySelector("#finding_capacitance").innerHTML = "C = " + capacitance + " F";
}

// Energy in Capacitor Calculation
function calc_energy_capacitor() {
    const C = parseFloat(document.getElementById('capacitance_energy').value);
    const V = parseFloat(document.getElementById('voltage_energy').value);

    if (isNaN(C) && isNaN(V)) {
        alert("Please enter valid values for Capacitance and Voltage.");
    }

    const capacitanceEnergy = Module.ccall(
        'energy_in_capacitor', 
        'number', 
        ['number', 'number'], 
        [C, V]
    );

    alert("U =  0.5 x C x V^2 \n U = 0.5 x " + C + "x" + V + "^2 \n U = " + capacitanceEnergy + " J");
    document.querySelector("#finding_energy").innerHTML = "U = " + capacitanceEnergy + " J";
}


// Modren Physics Calculation
// Energy of Photon Calculation
function calc_energy_photon() {
    const f = parseFloat(document.getElementById('frequency_photon').value);
    // const lambda = parseFloat(document.getElementById('wavelength_photon').value);

    if (isNaN(f)) {
        alert("Please enter valid values for Frequency.");
    } 

    const photonEnergy = Module.ccall(
        'energy_photon', 
        'number', 
        ['number', 'number'], 
        [f, 0]
    );

    alert("E = h x f \n E = 6.625 x 10^-34 x " + f + " \n E = " + photonEnergy + " J");
    document.querySelector("#finding_energy_photon").innerHTML = "E = " + photonEnergy + " J";
}

// de Broglie Wavelength Calculation
function calc_de_broglie_wavelength() {
    const m = parseFloat(document.getElementById('mass_broglie').value);
    const v = parseFloat(document.getElementById('velocity_broglie').value);

    if (isNaN(m) || isNaN(v)) {
        alert("Please enter valid values for Mass and Velocity.");
    }

    const debroglieWaveLength = Module.ccall(
        'de_broglie_wavelength', 
        'number', 
        ['number', 'number'], 
        [m, v]
    );
        alert("λ = h / mv \n λ = 6.625 x 10^-34 / " + m + " x " + v + " \n λ = " + debroglieWaveLength + " m");
    document.querySelector("#finding_de_broglie_wavelength").innerHTML = "λ = " + debroglieWaveLength + " m";
}

// Time Dilation Calculation
function calc_time_dilation() {
    const t0 = parseFloat(document.getElementById('proper_time_dilation').value);
    const v = parseFloat(document.getElementById('velocity_dilation').value);

    if (!isNaN(t0) && !isNaN(v)) {
        const timeDilation = Module.ccall(
            'time_dilation', 
            'number', 
            ['number', 'number'], 
            [t0, v]
        );

        alert("t = t0 / sqrt(1 - v^2 / c^2) \n t = " + t0 + " / sqrt(1 - " + v + "^2 / c^2) \n t = " + timeDilation + " s");
        document.querySelector("#finding_time_dilation").innerHTML = "t = " + timeDilation + " sec";
    } 
    else {
        alert("Please enter valid values for Proper Time and Velocity.");
    }
} 

// Length Contraction Calculation
function calc_length_contraction() {
    const l0 = parseFloat(document.getElementById('proper_length').value);
    const v = parseFloat(document.getElementById('velocity_length').value);

    if (!isNaN(l0) && !isNaN(v)) {
        const lengthContraction = Module.ccall(
            'length_contraction', 
            'number', 
            ['number', 'number'], 
            [l0, v]
        );

        alert("L = L0 / sqrt(1 - v^2 / c^2) \n L = " + l0 + " / sqrt(1 - " + v + "^2 / c^2) \n L = " + lengthContraction + " m");
        document.querySelector("#finding_length_contraction").innerHTML = "L = " + lengthContraction + " meters";
    } 
    else {
        alert("Please enter valid values for Proper Length and Velocity.");
    }
}


// Compton Wavelength Shift Calculation
// function calc_compton_shift() {
//     const theta = parseFloat(document.getElementById('angle_theta').value);

//     if (!isNaN(theta)) {
//         const result = Module.ccall('compton_wavelength_shift', 'number', ['number'], [theta]);
//         alert(`Compton Wavelength Shift (Δλ) is: ${result} meters`);
//     } else {
//         alert("Please enter a valid value for Scattering Angle.");
//     }
// }

// Nuclear Decay Calculation

// function calc_nuclear_decay() {
//     const N0 = parseFloat(document.getElementById('initial_quantity').value);
//     const lambda = parseFloat(document.getElementById('decay_constant').value);
//     const t = parseFloat(document.getElementById('time_decay').value);

//     if (!isNaN(N0) && !isNaN(lambda) && !isNaN(t)) {
//         const result = Module.ccall('nuclear_decay', 'number', ['number', 'number', 'number'], [N0, lambda, t]);
//         alert(`Remaining Quantity (N) is: ${result}`);
//     } else {
//         alert("Please enter valid values for Initial Quantity, Decay Constant, and Time.");
//     }
// }
