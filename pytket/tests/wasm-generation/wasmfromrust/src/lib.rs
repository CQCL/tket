// src/lib.rs

#[no_mangle]
pub extern "C" fn init() {    
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

