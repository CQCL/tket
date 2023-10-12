// src/lib.rs

#[no_mangle]
fn init() {    
}

#[no_mangle]
pub extern "C" fn add_one(x: i32) -> i32 {
    x + 1
}


#[no_mangle]
pub extern "C" fn multi(x: i32, y: i32) -> i32 {
    x * y
}


#[no_mangle]
pub extern "C" fn add_two(x: i32) -> i32 {
    x + 2
}

#[no_mangle]
pub extern "C" fn add_something(x: i64) -> i64 {
    x + 11
}

#[no_mangle]
pub extern "C" fn add_something_32(x: i32, y: i32) -> i32 {
    x + y
}

#[no_mangle]
pub extern "C" fn add_eleven(x: i32) -> i32 {
    x + 11
}

#[no_mangle]
pub extern "C" fn no_return(x: i32) {
    let _y = x + 11;
}

#[no_mangle]
pub extern "C" fn no_parameters() -> i32 {
    11
}

#[no_mangle]
pub extern "C" fn new_function() -> i32 {
    13
}

#[no_mangle]
pub extern "C" fn mixed_up(limit: i32) -> i32 {
    let mut i = 0;

    while i < limit {
        i = i * 2;
    }
return i
}

#[no_mangle]
pub fn mixed_up_2(limit: i32, limit2: i32) -> i32 {
    let mut i = 0;

    while i < limit {
        i = i * 2;
    }

    while i < limit2 {
        i = i * 3;
    }
return i
}


#[no_mangle]
fn mixed_up_3(limit: i32, limit2: i32, limit3: i32) -> i32 {
    let mut i = 0;

    while i < limit {
        i = i * 2;
    }

    while i < limit2 {
        i = i * 3;
    }

    while i < limit3 {
        i = i * 4;
    }

    return i
}


#[no_mangle]
fn unse_internal(p: i32) -> i32 {
    let mut r = no_parameters();

    r = add_eleven(r);
    r = add_something_32(r, p);

    return r
}
